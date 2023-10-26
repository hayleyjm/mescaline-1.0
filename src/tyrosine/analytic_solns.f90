module analytic_solns
  !
  ! A module that contains analytic solutions for the metrics in init_metrics for various things
  !   that are output by Mescaline
  !
  use options, only:dp,c_double,pi,ainit,HLinit,tinitial,&
       box_size,clen
  use prints, only:print_error
  use derivatives, only:deriv1fourth
  use periodic, only:apply_periodic,apply_periodic_fourth
  use manipulations, only:get_metric_at_pos
  implicit none


contains

    ! --------------------------------------------
    !
    ! FLRW scale factor a(time) in the longitudinal gauge
    !        --> time here is conformal time, eta
    !
    ! --------------------------------------------
    real(c_double) function aflrw(time)
        real(c_double), intent(in) :: time

        aflrw = ainit * (time/tinitial)**2

    end function aflrw

  
    ! --------------------------------------------
    !
    ! a subroutine to set the FLRW background at any time
    !
    ! --------------------------------------------
    subroutine set_bg_parameters(time,at,rhot,asq,rhostar,hub,adot,addot)
      real(c_double), intent(in) :: time
      real(c_double), intent(out) :: at,rhot,asq,rhostar,hub,adot,addot

      real(c_double) :: rhoinitial,hubinit,rhostar_tmp
      character(len=100) :: message, loc
      loc = "set_bg_parameters"

      call get_init_rho(rhoinitial,hubinit,rhostar)
      !
      ! set up arbitrary time background values
      at   = ainit * (time / tinitial)**2
      rhot = rhoinitial * ainit**3 / at**3
      hub  = sqrt(8._dp * pi * rhot * at**2 / 3._dp) ! Friedmann eqns.

      asq     = at * at  ! a^2
      adot    = hub * at ! a'
      addot   = 4._dp * pi * rhostar / 3._dp ! a''
      !
      ! ^^ NOTE that addot is POSITIVE for conformal time alp = a
      !     -- the minus sign in this eqn. is for proper time only!! (check RGTC!)
      !     -- this must be due to extra terms in R_mu,nu due to dt(g00) /neq 0
      !

      !
      ! check FLRW density is actually conserved
      rhostar_tmp = rhot * at**3
      if (abs(rhostar/rhostar_tmp - 1._dp)>1.e-15) then
         write(message,"(a,ES15.8,ES15.8)") " Conserved density NOT conserved: ",rhostar,rhostar_tmp
         call print_error(message,1,loc)
      endif

    end subroutine set_bg_parameters



    ! --------------------------------------------
    !
    ! get the initial density, hubble, and rhostar
    !
    ! --------------------------------------------
    subroutine get_init_rho(rhoinitial,hubinit,rhostar)
      real(c_double), intent(out) :: rhoinitial,hubinit,rhostar

      hubinit    = HLinit / box_size
      rhoinitial = 3._dp * hubinit**2 / (8._dp * pi * ainit**2)
      rhostar    = rhoinitial * ainit**3 ! always the same

    end subroutine get_init_rho






  !---------------------------------------------
  !
  ! Conformal hubble parameter + derivatives
  !
  !---------------------------------------------


  ! H'
  real(c_double) function hubdash(hub,aval,addot)
    real(c_double) :: hub,addot,aval

    hubdash = addot / aval  - hub**2

  end function hubdash

  ! H''
  real(c_double) function hubddash(hub,aval,addot)
    real(c_double) :: hub,addot,aval

    !hubddash = hub * addot / aval - 2._dp*hub**3
    hubddash = 2._dp*hub**3 - 3._dp * hub * addot / aval

  end function hubddash

  !---------------------------------------------
  !---------------------------------------------

  !
  ! function to return the EXACT EdS dL(z)
  real(c_double) function EdS_dL(H0,z)
      real(c_double), intent(in) :: H0,z
      real(c_double) :: bracfac

      bracfac = 1._dp + z - sqrt(1._dp + z)
      EdS_dL = 2._dp*bracfac/H0

  end function EdS_dL


  !
  ! function to return the LCDM/FLRW DL(z) accurate to third order in redshift
  !      for FLAT only -- forced to EdS for now
  real(c_double) function LCDM_dL3(H0,z)
      real(c_double), intent(in) :: z,H0 ! redshift, Hubble(z=0) = H0
      real(c_double) :: q0,dL1,dL2,dL3,oL,om
      ! Force to EdS
      oL = 0._dp   ! omega_lambda
      om = 1.0_dp  ! omega_m

      q0  = 0.5_dp*om - oL ! deceleration parameter q_0 (https://people.ast.cam.ac.uk/~pettini/Intro%20Cosmology/Lecture06.pdf)
      dL1 = 1._dp/H0
      dL2 = 0.5_dp*(1._dp - q0)/H0
      dL3 = (3._dp*q0**2 + q0 - 2._dp)/(6._dp*H0)

      LCDM_dL3 = dL1*z + dL2*z**2 + dL3*z**3

  end function LCDM_dL3


  ! --------------------------------------------
  !
  ! Conformal scaled time \xi
  !
  ! --------------------------------------------
  real(c_double) function xi(time,rhostar)
    real(c_double) :: time,rhostar

    xi = 1._dp + sqrt(2._dp * pi * rhostar / (3._dp * ainit) ) * (time - tinitial)

  end function xi


  !
  ! a subroutine to get the time at z=0 given some data set in options
  !
  ! this comes from z = (1+zini)/a - 1
  !      then sub in a=aini(t/tini)^2 and rearrange setting z=0
  !
  real(c_double) function etaz0(zinit)
      real(c_double) :: zinit

      etaz0 = tinitial * sqrt((1._dp + zinit) / ainit)

  end function etaz0


  ! --------------------------------------------
  !
  ! Exact solutions for the FLRW metric
  !
  ! --------------------------------------------
  subroutine FLRW_analytic(nx,time,aD,hub,traceR,traceK,Theta,sigxx,sigma2,w2)
    integer, intent(in) :: nx
    real(c_double), intent(in) :: time
    real(c_double), intent(out) :: aD,hub
    real(c_double), intent(out), dimension(nx,nx,nx) :: traceR,traceK,Theta,sigma2,w2,sigxx ! one component of sigmadd

    real(c_double) :: rhot,asq,rhostar,adot,addot

    !
    ! some things are zero for FLRW
    traceR = 0._dp; sigma2 = 0._dp
    sigxx = 0._dp; w2 = 0._dp

    !
    ! some things are not
    call set_bg_parameters(time,aD,rhot,asq,rhostar,hub,adot,addot)

    Theta  = 3._dp * adot / asq
    traceK = - Theta

  end subroutine FLRW_analytic




end module analytic_solns
