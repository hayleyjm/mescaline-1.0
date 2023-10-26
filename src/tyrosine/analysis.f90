!-------------------------------------------------------------------------
!
!  Module that sends data for the actual post-simulation analysis
!
!-------------------------------------------------------------------------
module analysis
  use options,   only:c_double,clen,domain_type
  use prints,    only:print_info,print_error
  use cactoidae, only: gotalp,gotrho
  implicit none

contains

  !
  ! This subroutine takes data at ntimes and send for analysis (tiny little baby routine)
  !
  subroutine send_for_analysis(nx,ntimes,times,xmax,xmin,dx,gij,kij,alp,vel0,vel1,vel2,rho,dens)
    use options, only:nspheres,radius,looprad,radmin,radmax,drad
    use ricci,   only:compute_ricci

    integer, intent(in) :: nx,ntimes
    real(c_double), dimension(6,nx,nx,nx,ntimes),  intent(in) :: gij,kij
    real(c_double), dimension(nx,nx,nx,ntimes),    intent(in) :: vel0,vel1,vel2
    real(c_double), dimension(nx,nx,nx),    intent(in) :: dens
    real(c_double), dimension(nx,nx,nx), intent(inout) :: rho ! inout because we need this for ricci
    real(c_double), intent(in) :: times(ntimes),xmax(3),xmin(3),dx,alp(nx,nx,nx)

    real(c_double) :: rad,radend
    integer :: nsph
    character(len=clen) :: message,loc
    loc = " send_for_analysis" ! our current location for print messages

    !
    ! First some information about the domain type to be clear
    !
    select case(domain_type)
    case("sphere")
        call print_info("Averaging over spherical domains.",loc)
    case("cube")
        call print_info("Averaging over cubic domains (I might refer to them as 'spheres' somewhere, sorry, &
        & I'm working on this. They are cubes though, I promise)",loc)
    case default
        call print_error("Please choose a valid domain_type",2,loc)
    end select

    if (looprad) then
       write(message,"(a,f7.1,a,f7.1,a,f7.1)") " Looping from r_D = ",radmin," to r_D = ",radmax," with dr_D = ",drad
       call print_info(message, loc)
       rad    = radmin
       radend = radmax
    else
       rad    = radius ! single radius of averaging
       radend = radius
    endif

    do while (rad<=radend)
       nsph = nspheres
       call check_radius(xmax,xmin,dx,rad,nsph)
       call compute_ricci(nx,ntimes,dx,xmax(1),xmin(1),rad,nsph,gij,kij,rho,&
            & vel0,vel1,vel2,alp,times,gotrho,gotalp,dens)
       if (looprad) then
          rad = rad + drad
       else
          exit
       endif
    enddo

  end subroutine send_for_analysis



  !
  ! check radius of averaging for given box size
  !
  subroutine check_radius(xmax,xmin,dx,rad,nspheres)
    real(c_double), intent(in) :: xmax(3), xmin(3), dx
    real(c_double), intent(inout) :: rad
    integer, intent(inout) :: nspheres
    real(c_double) :: halfbox
    character(len=clen) :: message,loc
    loc = "check_radius"

    halfbox = (xmax(1) - xmin(1)) / 2.d0
    !
    ! Check our chosen radius relative to grid spacing
    !
    if (rad/=0.) then
       !
       ! check these if we are NOT all-box
       !
       if (rad<=dx) then
          call print_error( "Domain radius < grid spacing. This isn't going to work. ",2,loc)
       elseif (rad>dx .and. rad<=2.*dx) then
          call print_error( "Domain radius < 2 dx. Please increase this. ",1,loc)
       elseif (rad>dx .and. rad<=5.*dx) then
          call print_error( "Domain radius < 5 dx. Not ideal. ",1,loc)
       endif

    endif

    !
    ! Check our chosen radius vs. number of spheres
    !
    if (rad==0. .and. nspheres/=1) then
       !
       ! We chose an all-box average but more than one sphere - makes no sense silly
       !
       call print_error(" OVERRIDING nspheres = 1 since rad = 0 ",0,loc)
       nspheres = 1
    endif
    if (rad==halfbox .and. nspheres/=1) then
       !
       ! We have sphere/s with radius the size of the box
       ! Set nspheres=1 to be sure
       !
       call print_error(" OVERRIDING nspheres = 1 since radius is half the box size ",0,loc)
       nspheres = 1
    elseif (rad>halfbox) then
       !
       ! You set the sphere larger than the box -- overriding
       !
       nspheres = 1
       write(message,"(a,f7.1)") " You can't have radius = ",rad
       call print_error(message,0,loc)
       write(message,"(a,f7.1,a,i6)") "OVERRIDING radius = ",&
        & halfbox," (boxL/2) and nspheres = ",nspheres
       call print_error(message,0,loc)
       rad = halfbox
    end if

  end subroutine check_radius


end module analysis
