#include "fabm_driver.h"

module pisces_pom_remineralization
   use fabm_types
   use fabm_particle
   use pisces_shared

   implicit none

   private

   type, extends(type_particle_model), public :: type_pisces_pom_remineralization
      type (type_state_variable_id) :: id_c, id_fe, id_poc, id_sfe, id_doc, id_fer
      type (type_dependency_id)     :: id_zremi, id_ws
      real(rk) :: solgoc
   contains
      procedure :: initialize => pom_remineralization_initialize
      procedure :: do         => pom_remineralization_do
   end type

   type, extends(type_base_model) :: type_pisces_lability
      type (type_surface_dependency_id) :: id_hmld
      type (type_dependency_id) :: id_tem, id_gdept_n, id_e3t_n, id_ws, id_cons, id_prod, id_poc
      type (type_diagnostic_variable_id) :: id_zremi
      logical :: mlvar
      integer :: jcpoc
      real(rk) :: xremip, rshape
      real(rk), allocatable :: alphan(:), reminp(:)
   contains
      procedure :: initialize => lability_initialize
      procedure :: do_column  => lability_do_column
   end type

contains

   subroutine pom_remineralization_initialize(self, configunit)
      class (type_pisces_pom_remineralization), intent(inout), target :: self
      integer,                                  intent(in)            :: configunit

      class (type_pisces_lability), pointer :: lability

      allocate(lability)

      ! Here we compute the GOC -> POC rate due to the shrinking
      ! of the fecal pellets/aggregates as a result of bacterial
      ! solubilization
      ! This is based on a fractal dimension of 2.56 and a spectral
      ! slope of -3.6 (identical to what is used in p4zsink to compute
      ! aggregation
      call self%get_parameter(self%solgoc, 'solgoc', '-', 'fraction broken down to smaller POM', default=0.04_rk/ 2.56 * 1./ ( 1.-50**(-0.04) ))
      call self%get_parameter(lability%mlvar, 'mlvar', '-', 'lability varies in mixing layer', default=.false.)
      call self%get_parameter(lability%xremip, 'xremip', '-', 'remineralisation rate', default=0.035_rk)
      call self%get_parameter(lability%jcpoc, 'jcpoc', '', 'number of lability classes', default=15)
      call self%get_parameter(lability%rshape, 'rshape', '-', 'shape of the gamma function', default=1._rk)

      call self%add_child(lability, 'lability', configunit=-1)
      call lability%request_coupling(lability%id_poc, '../c')
      call lability%request_coupling(lability%id_ws, '../ws')

      call self%register_state_dependency(self%id_c, 'c', 'mol C L-1', 'particulate organic carbon')
      call self%register_state_dependency(self%id_fe, 'fe', 'mol C L-1', 'particulate organic carbon')
      call self%register_state_dependency(self%id_doc, 'doc', 'mol C L-1', 'dissolved organic carbon')
      call self%register_state_dependency(self%id_fer, 'fer', 'mol C L-1', 'iron')
      call self%register_state_dependency(self%id_poc, 'poc', 'mol C L-1', 'small particulate organic carbon')
      call self%register_state_dependency(self%id_sfe, 'sfe', 'mol C L-1', 'small particulate organic iron')
      call self%register_dependency(self%id_ws, 'ws', 'm d-1', 'sinking velocity')
      call self%register_dependency(self%id_zremi, 'zremi', 's-1', 'remineralization rate')

      call self%request_coupling_to_model(self%id_c, 'pom', 'c')
      call self%request_coupling_to_model(self%id_fe, 'pom', 'fe')
      call self%request_coupling_to_model(self%id_ws, 'pom', 'ws')
      call self%request_coupling_to_model(self%id_doc, 'dom', 'c')
      if (self%solgoc == 0._rk) then
         call self%request_coupling(self%id_poc, 'zero')
         call self%request_coupling(self%id_sfe, 'zero')
      else
         call self%request_coupling_to_model(self%id_poc, 'spom', 'c')
         call self%request_coupling_to_model(self%id_sfe, 'spom', 'fe')
      end if
      call self%request_coupling(self%id_zremi, 'lability/zremi')
   end subroutine

   subroutine pom_remineralization_do(self, _ARGUMENTS_DO_)
      class (type_pisces_pom_remineralization), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: c, fe, zremi, zorem, zofer

      _LOOP_BEGIN_
         _GET_(self%id_c, c)
         _GET_(self%id_fe, fe)
         _GET_(self%id_zremi, zremi)
         zorem = zremi * c
         zofer = zremi * fe
         _ADD_SOURCE_(self%id_c,   - (1._rk + self%solgoc) * zorem)
         _ADD_SOURCE_(self%id_fe,  - (1._rk + self%solgoc) * zofer)
         _ADD_SOURCE_(self%id_doc, + zorem)
         _ADD_SOURCE_(self%id_fer, + zofer)
         _ADD_SOURCE_(self%id_poc, + self%solgoc * zorem)
         _ADD_SOURCE_(self%id_sfe, + self%solgoc * zofer)
      _LOOP_END_
   end subroutine

   subroutine lability_initialize(self, configunit)
      class (type_pisces_lability), intent(inout), target :: self
      integer,                      intent(in)            :: configunit

      integer :: jn, ifault
      real(rk) :: remindelta, reminup, remindown

      call self%register_diagnostic_variable(self%id_zremi, 'zremi', 's-1', 'remineralization rate', source=source_do_column)

      call self%register_dependency(self%id_ws, 'ws', 'm d-1', 'sinking velocity of particulate organic silicon')
      call self%register_dependency(self%id_e3t_n, standard_variables%cell_thickness)
      call self%register_dependency(self%id_gdept_n, standard_variables%depth)
      call self%register_dependency(self%id_hmld, turbocline_depth)
      call self%register_dependency(self%id_tem, standard_variables%temperature)
      call self%register_dependency(self%id_poc, 'poc', 'mol C L-1', 'particulate organic carbon')
      call self%register_dependency(self%id_prod, 'prod', 'mol C L-1 s-1', 'particulate organic carbon sources')
      call self%register_dependency(self%id_cons, 'cons', 'mol C L-1 s-1', 'particulate organic carbon sinks')
      call self%request_coupling(self%id_prod, 'zero')
      call self%request_coupling(self%id_cons, 'zero')

      ALLOCATE( self%alphan(self%jcpoc) , self%reminp(self%jcpoc) )
      !
      IF (self%jcpoc > 1) THEN
         !
         remindelta = LOG(4. * 1000. ) / REAL(self%jcpoc-1, rk)
         reminup = 1./ 400. * EXP(remindelta)
         !
         ! Discretization based on incomplete gamma functions
         ! As incomplete gamma functions are not available in standard 
         ! fortran 95, they have been coded as functions in this module (gamain)
         ! ---------------------------------------------------------------------
         !
         self%alphan(1) = gamain(reminup, self%rshape, ifault)
         self%reminp(1) = gamain(reminup, self%rshape+1.0, ifault) * self%xremip / self%alphan(1)
         DO jn = 2, self%jcpoc-1
            reminup = 1./ 400. * EXP( REAL(jn, rk) * remindelta)
            remindown = 1. / 400. * EXP( REAL(jn-1, rk) * remindelta)
            self%alphan(jn) = gamain(reminup, self%rshape, ifault) - gamain(remindown, self%rshape, ifault)
            self%reminp(jn) = gamain(reminup, self%rshape+1.0, ifault) - gamain(remindown, self%rshape+1.0, ifault)
            self%reminp(jn) = self%reminp(jn) * self%xremip / self%alphan(jn)
         END DO
         remindown = 1. / 400. * EXP( REAL(self%jcpoc-1, rk) * remindelta)
         self%alphan(self%jcpoc) = 1.0 - gamain(remindown, self%rshape, ifault)
         self%reminp(self%jcpoc) = 1.0 - gamain(remindown, self%rshape+1.0, ifault)
         self%reminp(self%jcpoc) = self%reminp(self%jcpoc) * self%xremip / self%alphan(self%jcpoc)

      ELSE
         self%alphan(self%jcpoc) = 1.
         self%reminp(self%jcpoc) = self%xremip
      ENDIF
   end subroutine

   subroutine lability_do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_pisces_lability), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: poc, zremi, zdep, gdept_n, e3t_n, ws, tem, cons, prod
      real(rk) :: totprod, totthick, totcons
      real(rk) :: e3t_n_prev, gdept_n_prev, ws_prev, tgfunc_prev, poc_prev, cons_prev, prod_prev
      logical :: first
      integer :: jn
      real(rk) :: tgfunc, alphat, remint, zsizek, zsizek1, zpoc, ztremint
      real(rk) :: alpha(self%jcpoc), alpha_prev(self%jcpoc)

      _GET_SURFACE_(self%id_hmld, zdep)

      ztremint = self%xremip
      alpha = self%alphan

      IF (self%mlvar) THEN
         totprod = 0._rk
         totthick = 0._rk
         totcons = 0._rk
         _DOWNWARD_LOOP_BEGIN_
            _GET_(self%id_gdept_n, gdept_n)
            IF (gdept_n <= zdep ) THEN
               _GET_(self%id_e3t_n, e3t_n)
               _GET_(self%id_cons, cons)
               _GET_(self%id_prod, prod)
               _GET_(self%id_tem, tem)
               _GET_(self%id_poc, poc)
               _GET_(self%id_ws, ws)
               tgfunc = EXP( 0.063913_rk * tem )  ! Jorn: Eq 4a in PISCES-v2 paper, NB EXP(0.063913) = 1.066 = b_P
               totprod = totprod + prod * e3t_n * rday
               totthick = totthick + e3t_n * tgfunc
               totcons = totcons - cons * e3t_n * rday / ( poc + rtrn )
            ENDIF
         _DOWNWARD_LOOP_END_

         alphat = 0.0
         remint = 0.0
         DO jn = 1, self%jcpoc
            ! For each lability class, the system is supposed to be 
            ! at equilibrium: Prod - Sink - w alphap = 0.
            alpha(jn) = totprod * self%alphan(jn) / ( self%reminp(jn)    &
            &                     * totthick + totcons + ws + rtrn )
            alphat = alphat + alpha(jn)
         END DO
         DO jn = 1, self%jcpoc
            alpha(jn) = alpha(jn) / ( alphat + rtrn)
            remint = remint + alpha(jn) * self%reminp(jn)
         END DO
         ! Mean remineralization rate in the mixed layer
         ztremint =  MAX( 0., remint )
      ENDIF

      first = .true.

      _DOWNWARD_LOOP_BEGIN_
         _GET_(self%id_gdept_n, gdept_n)
         _GET_(self%id_e3t_n, e3t_n)
         _GET_(self%id_tem, tem)
         _GET_(self%id_ws, ws)
         _GET_(self%id_poc, poc)
         _GET_(self%id_cons, cons)
         _GET_(self%id_prod, prod)

         tgfunc = EXP( 0.063913_rk * tem )  ! Jorn: Eq 4a in PISCES-v2 paper, NB EXP(0.063913) = 1.066 = b_P
         zsizek = e3t_n / 2. / (ws + rtrn) * tgfunc

         !
         ! In the case of GOC, lability is constant in the mixed layer 
         ! It is computed only below the mixed layer depth
         ! ------------------------------------------------------------
         !
         IF( gdept_n > zdep .and. .not. first) THEN
            alphat = 0.
            remint = 0.
            !
            !
            IF ( gdept_n_prev <= zdep ) THEN
               ! 
               ! The first level just below the mixed layer needs a 
               ! specific treatment because lability is supposed constant
               ! everywhere within the mixed layer. This means that 
               ! change in lability in the bottom part of the previous cell
               ! should not be computed
               ! ----------------------------------------------------------
               !
               ! POC concentration is computed using the lagrangian 
               ! framework. It is only used for the lability param
               zpoc = poc_prev + cons * rday               &   ! Jorn: dropped division by rfact2 as consgoc is already per second (unlike in pisces, where it is premultiplied by time step rfact2)
               &   * e3t_n / 2. / (ws + rtrn)                     ! Jorn: time (d) it takes to sink half the current layer height?
               zpoc = MAX(0., zpoc)
               !
               DO jn = 1, self%jcpoc
                  !
                  ! Lagrangian based algorithm. The fraction of each 
                  ! lability class is computed starting from the previous
                  ! level
                  ! -----------------------------------------------------
                  !
                  ! the concentration of each lability class is calculated
                  ! as the sum of the different sources and sinks
                  ! Please note that production of new GOC experiences
                  ! degradation 
                  alpha(jn) = alpha_prev(jn) * exp( -self%reminp(jn) * zsizek ) * zpoc &
                  &   + prod * self%alphan(jn) / tgfunc / self%reminp(jn)             &  ! Jorn: same for POC, except that GOC->POC conversion [zorem3*alphag] is added to prodgoc*alphan
                  &   * ( 1. - exp( -self%reminp(jn) * zsizek ) ) * rday   ! Jorn: dropped division by rfact2 as prodgoc is already per second (unlike in pisces, where it is premultiplied by time step rfact2)
                  alphat = alphat + alpha(jn)
                  remint = remint + alpha(jn) * self%reminp(jn)
               END DO
            ELSE
               !
               ! standard algorithm in the rest of the water column
               ! See the comments in the previous block.
               ! ---------------------------------------------------
               !
               zpoc = poc_prev + cons_prev * rday               &   ! Jorn: dropped division by rfact2 as consgoc is already per second (unlike in pisces, where it is premultiplied by time step rfact2)
               &   * e3t_n_prev / 2. / (ws_prev + rtrn) + cons   &
               &   * rday * e3t_n / 2. / (ws + rtrn)
               zpoc = max(0., zpoc)
               !
               DO jn = 1, self%jcpoc
                  alpha(jn) = alpha_prev(jn) * exp( -self%reminp(jn) * ( zsizek              &
                  &   + zsizek1 ) ) * zpoc + ( prod_prev / tgfunc_prev * ( 1.           &
                  &   - exp( -self%reminp(jn) * zsizek1 ) ) * exp( -self%reminp(jn) * zsizek ) + prod &  ! Jorn: same for POC, except that GOC->POC conversion [zorem3*alphag] is added to prodgoc*alphan
                  &   / tgfunc * ( 1. - exp( -self%reminp(jn) * zsizek ) ) ) * rday / self%reminp(jn) * self%alphan(jn)    ! Jorn: dropped division by rfact2 as prodgoc is already per second (unlike in pisces, where it is premultiplied by time step rfact2)
                  alphat = alphat + alpha(jn)
                  remint = remint + alpha(jn) * self%reminp(jn)
               END DO
            ENDIF
            !
            DO jn = 1, self%jcpoc
               ! The contribution of each lability class at the current
               ! level is computed
               alpha(jn) = alpha(jn) / ( alphat + rtrn)
            END DO
            ! Computation of the mean remineralisation rate
            ztremint =  MAX(0., remint / ( alphat + rtrn) )
         ENDIF
         zremi = MIN( self%xremip , ztremint )
         _SET_DIAGNOSTIC_(self%id_zremi, zremi * xstep * tgfunc)

         e3t_n_prev = e3t_n
         gdept_n_prev = gdept_n
         ws_prev = ws
         tgfunc_prev = tgfunc
         poc_prev = poc
         zsizek1  = zsizek
         alpha_prev = alpha
         cons_prev = cons
         prod_prev = prod
         first = .false.
      _DOWNWARD_LOOP_END_
   end subroutine

   REAL FUNCTION alngam( xvalue, ifault )
      !*****************************************************************************80
      !
      !! ALNGAM computes the logarithm of the gamma function.
      !
      !  Modified:    13 January 2008
      !
      !  Author  :    Allan Macleod
      !               FORTRAN90 version by John Burkardt
      !
      !  Reference:
      !    Allan Macleod, Algorithm AS 245,
      !    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
      !    Applied Statistics,
      !    Volume 38, Number 2, 1989, pages 397-402.
      !
      !  Parameters:
      !
      !    Input, real ( kind = 8 ) XVALUE, the argument of the Gamma function.
      !
      !    Output, integer ( kind = 4 ) IFAULT, error flag.
      !    0, no error occurred.
      !    1, XVALUE is less than or equal to 0.
      !    2, XVALUE is too big.
      !
      !    Output, real ( kind = 8 ) ALNGAM, the logarithm of the gamma function of X.
      !*****************************************************************************80
  implicit none

  real(rk), parameter :: alr2pi = 0.918938533204673E+00
  integer:: ifault
  real(rk), dimension ( 9 ) :: r1 = (/ &
    -2.66685511495E+00, &
    -24.4387534237E+00, &
    -21.9698958928E+00, &
     11.1667541262E+00, &
     3.13060547623E+00, &
     0.607771387771E+00, &
     11.9400905721E+00, &
     31.4690115749E+00, &
     15.2346874070E+00 /)
  real(rk), dimension ( 9 ) :: r2 = (/ &
    -78.3359299449E+00, &
    -142.046296688E+00, &
     137.519416416E+00, &
     78.6994924154E+00, &
     4.16438922228E+00, &
     47.0668766060E+00, &
     313.399215894E+00, &
     263.505074721E+00, &
     43.3400022514E+00 /)
  real(rk), dimension ( 9 ) :: r3 = (/ &
    -2.12159572323E+05, &
     2.30661510616E+05, &
     2.74647644705E+04, &
    -4.02621119975E+04, &
    -2.29660729780E+03, &
    -1.16328495004E+05, &
    -1.46025937511E+05, &
    -2.42357409629E+04, &
    -5.70691009324E+02 /)
  real(rk), dimension ( 5 ) :: r4 = (/ &
     0.279195317918525E+00, &
     0.4917317610505968E+00, &
     0.0692910599291889E+00, &
     3.350343815022304E+00, &
     6.012459259764103E+00 /)
  real (rk) :: x
  real (rk) :: x1
  real (rk) :: x2
  real (rk), parameter :: xlge = 5.10E+05
  real (rk), parameter :: xlgst = 1.0E+30
  real (rk) :: xvalue
  real (rk) :: y

  x = xvalue
  alngam = 0.0E+00
!
!  Check the input.
!
  if ( xlgst <= x ) then
    ifault = 2
    return
  end if
  if ( x <= 0.0E+00 ) then
    ifault = 1
    return
  end if

  ifault = 0
!
!  Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
!
  if ( x < 1.5E+00 ) then

    if ( x < 0.5E+00 ) then
      alngam = - log ( x )
      y = x + 1.0E+00
!
!  Test whether X < machine epsilon.
!
      if ( y == 1.0E+00 ) then
        return
      end if

    else

      alngam = 0.0E+00
      y = x
      x = ( x - 0.5E+00 ) - 0.5E+00

    end if

    alngam = alngam + x * (((( &
        r1(5)   * y &
      + r1(4) ) * y &
      + r1(3) ) * y &
      + r1(2) ) * y &
      + r1(1) ) / (((( &
                  y &
      + r1(9) ) * y &
      + r1(8) ) * y &
      + r1(7) ) * y &
      + r1(6) )

    return

  end if
!
!  Calculation for 1.5 <= X < 4.0.
!
  if ( x < 4.0E+00 ) then

    y = ( x - 1.0E+00 ) - 1.0E+00

    alngam = y * (((( &
        r2(5)   * x &
      + r2(4) ) * x &
      + r2(3) ) * x &
      + r2(2) ) * x &
      + r2(1) ) / (((( &
                  x &
      + r2(9) ) * x &
      + r2(8) ) * x &
      + r2(7) ) * x &
      + r2(6) )
!
!  Calculation for 4.0 <= X < 12.0.
!
  else if ( x < 12.0E+00 ) then

    alngam = (((( &
        r3(5)   * x &
      + r3(4) ) * x &
      + r3(3) ) * x &
      + r3(2) ) * x &
      + r3(1) ) / (((( &
                  x &
      + r3(9) ) * x &
      + r3(8) ) * x &
      + r3(7) ) * x &
      + r3(6) )
!
!  Calculation for 12.0 <= X.
!
  else

    y = log ( x )
    alngam = x * ( y - 1.0E+00 ) - 0.5E+00 * y + alr2pi

    if ( x <= xlge ) then

      x1 = 1.0E+00 / x
      x2 = x1 * x1

      alngam = alngam + x1 * ( ( &
             r4(3)   * &
        x2 + r4(2) ) * &
        x2 + r4(1) ) / ( ( &
        x2 + r4(5) ) * &
        x2 + r4(4) )

    end if

   end if

   END FUNCTION alngam


   REAL FUNCTION gamain( x, p, ifault )
!*****************************************************************************80
!
!! GAMAIN computes the incomplete gamma ratio.
!
!  Discussion:
!
!    A series expansion is used if P > X or X <= 1.  Otherwise, a
!    continued fraction approximation is used.
!
!  Modified:
!
!    17 January 2008
!
!  Author:
!
!    G Bhattacharjee
!    FORTRAN90 version by John Burkardt
!
!  Reference:
!
!    G Bhattacharjee,
!    Algorithm AS 32:
!    The Incomplete Gamma Integral,
!    Applied Statistics,
!    Volume 19, Number 3, 1970, pages 285-287.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, P, the parameters of the incomplete 
!    gamma ratio.  0 <= X, and 0 < P.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no errors.
!    1, P <= 0.
!    2, X < 0.
!    3, underflow.
!    4, error return from the Log Gamma routine.
!
!    Output, real ( kind = 8 ) GAMAIN, the value of the incomplete
!    gamma ratio.
!
  implicit none

  real (rk) a
  real (rk), parameter :: acu = 1.0E-08
  real (rk) an
  real (rk) arg
  real (rk) b
  real (rk) dif
  real (rk) factor
  real (rk) g
  real (rk) gin
  integer i
  integer ifault
  real (rk), parameter :: oflo = 1.0E+37
  real (rk) p
  real (rk) pn(6)
  real (rk) rn
  real (rk) term
  real (rk), parameter :: uflo = 1.0E-37
  real (rk) x
!
!  Check the input.
!
  if ( p <= 0.0E+00 ) then
    ifault = 1
    gamain = 0.0E+00
    return
  end if

  if ( x < 0.0E+00 ) then
    ifault = 2
    gamain = 0.0E+00
    return
  end if

  if ( x == 0.0E+00 ) then
    ifault = 0
    gamain = 0.0E+00
    return
  end if

  g = alngam ( p, ifault )

  if ( ifault /= 0 ) then
    ifault = 4
    gamain = 0.0E+00
    return
  end if

  arg = p * log ( x ) - x - g

  if ( arg < log ( uflo ) ) then
    ifault = 3
    gamain = 0.0E+00
    return
  end if

  ifault = 0
  factor = exp ( arg )
!
!  Calculation by series expansion.
!
  if ( x <= 1.0E+00 .or. x < p ) then

    gin = 1.0E+00
    term = 1.0E+00
    rn = p

    do

      rn = rn + 1.0E+00
      term = term * x / rn
      gin = gin + term

      if ( term <= acu ) then
        exit
      end if

    end do

    gamain = gin * factor / p
    return

  end if
!
!  Calculation by continued fraction.
!
  a = 1.0E+00 - p
  b = a + x + 1.0E+00
  term = 0.0E+00

  pn(1) = 1.0E+00
  pn(2) = x
  pn(3) = x + 1.0E+00
  pn(4) = x * b

  gin = pn(3) / pn(4)

  do

    a = a + 1.0E+00
    b = b + 2.0E+00
    term = term + 1.0E+00
    an = a * term
    do i = 1, 2
      pn(i+4) = b * pn(i+2) - an * pn(i)
    end do

    if ( pn(6) /= 0.0E+00 ) then

      rn = pn(5) / pn(6)
      dif = abs ( gin - rn )
!
!  Absolute error tolerance satisfied?
!
      if ( dif <= acu ) then
!
!  Relative error tolerance satisfied?
!
        if ( dif <= acu * rn ) then
          gamain = 1.0E+00 - factor * gin
          exit
        end if

      end if

      gin = rn

    end if

    do i = 1, 4
      pn(i) = pn(i+2)
    end do
    if ( oflo <= abs ( pn(5) ) ) then

      do i = 1, 4
        pn(i) = pn(i) / oflo
      end do

    end if

  end do

   END FUNCTION gamain

end module