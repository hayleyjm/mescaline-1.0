module spines
  !
  ! A module containing the core calculations for spatial averaging
  !
  use options, only:pi,c_double,domain_type,clen
  use prints, only:print_error
  implicit none

contains

  ! ======================================================================
  !
  !                        -------- AVERAGING ---------
  !
  ! ======================================================================

  !
  ! Calculate averages over a given radius for *one sphere* only (OR cube)
  !
  subroutine calc_average_within_radius(xi,xj,xk,origin,radius,detg,dx,avg,posval)
    real(c_double), intent(inout) :: avg
    real(c_double), intent(in) :: radius, dx, detg
    real(c_double), intent(in) :: posval  ! quantity to be averaged at position (xi,xj,xk)
    real(c_double), intent(in) :: xi,xj,xk
    real(c_double), intent(in) :: origin(3)
    real(c_double) :: rval2,xrel,yrel,zrel
    character(len=clen) :: loc
    loc = " calc_average_within_radius"

    if (radius==0.) then
       !
       ! Calculate average over whole box
       !
       avg = avg + sqrt(detg) * posval * dx**3
    else
       !
       ! Calculate where we are in radius
       !
       ! 1. Get coordinates relative to origin of domain
       xrel = abs(xi - origin(1))
       yrel = abs(xj - origin(2))
       zrel = abs(xk - origin(3))
       ! 2. Check if we are inside this domain
       select case(domain_type)
       case("sphere")
           !
           ! Spherical domain; check if we are within radius
           rval2 = xrel**2 + yrel**2 + zrel**2
           if (rval2<=radius**2) then
               !
               ! We are within the radius of averaging, add to the average
               !
               avg = avg + sqrt(detg) * posval * dx**3
           endif
       case("cube")
           !
           ! Cubic domain; check that ALL x,y,z are < rad
           if (xrel<=radius .and. yrel<=radius .and. zrel<=radius) then
               !
               ! We are within the cubic domain, add to the average
               !
               avg = avg + sqrt(detg) * posval * dx**3
           endif

       case default
           call print_error("Please choose a valid domain_type",2,loc)
       end select

    endif

  end subroutine calc_average_within_radius


  !
  ! Calculate averages within a given radius, for MANY spheres!
  !        --> calculates avgs over all the variables we want at once
  !
  subroutine calc_average_within_radius_manyspheres(xi,xj,xk,radius,nspheres,randorigins,detg,dx,navgs,avgs,posvals)
    integer, intent(in) :: nspheres, navgs        ! number of spheres, number of variables we're averaging
    real(c_double), intent(in) :: xi,xj,xk        ! current position in x,y,z-values
    real(c_double), intent(in) :: radius, dx, detg
    real(c_double), intent(inout) :: avgs(navgs,nspheres)  ! list of average arrays of size=numspheres
    real(c_double), intent(in) :: posvals(navgs)           ! values to add to average at current position
    real(c_double), intent(in) :: randorigins(3,nspheres)  ! origins of the numspheres spheres we have

    real(c_double) :: rval2,radius2,origin(3),sdetg,xrel,yrel,zrel,radiusp
    real(c_double), parameter :: abit = 1.e-5 ! a bit to add to radius to make sure we get the whole thing
    integer :: l,k
    character(len=clen) :: loc
    loc = " calc_average_within_radius_manyspheres"

    radiusp = radius + abit
    radius2 = radiusp*radiusp
    sdetg   = sqrt(detg)

    !print*,navgs
    !avgs(3,:),posvals(3)

    if (radius==0.) then
       !
       ! Calculate avg over whole box
       !
       do k=1,navgs ! Loop over variables
          avgs(k,:) = avgs(k,:) + sdetg * dx**3 * posvals(k)
       enddo
    else
       !
       ! Calculate where we are in radius for each domain position
       !
       do l=1,nspheres
          origin = randorigins(:,l)
          ! 1. Get distance from origin of domain
          xrel = abs(xi - origin(1))
          yrel = abs(xj - origin(2))
          zrel = abs(xk - origin(3))
          ! 2. Check if we are inside this domain
          select case(domain_type)
          case("sphere")
              !
              ! Spherical domain; check if we are within radius
              rval2 = xrel**2 + yrel**2 + zrel**2
              if (rval2<=radius2) then
                  !
                  ! We are within this sphere, add to average
                  !
                  do k=1,navgs ! Loop over variables
                      avgs(k,l) = avgs(k,l) + sdetg * dx**3 * posvals(k)
                  enddo
              endif
              !
              ! Else: not in this sphere, move on and check the next one
              !
          case("cube")
              !
              ! Cubic domain; check that ALL x,y,z are < rad
              if (xrel<=radiusp .and. yrel<=radiusp .and. zrel<=radiusp) then
                  !
                  ! We are within the cubic domain, add to the average
                  !
                  do k=1,navgs ! Loop over variables
                      avgs(k,l) = avgs(k,l) + sdetg * dx**3 * posvals(k)
                  enddo
              endif

          case default
              call print_error("Please choose a valid domain_type",2,loc)
          end select
       enddo
    endif

  end subroutine calc_average_within_radius_manyspheres


  !
  ! Helps calculate the fraction of OD vs UD regions
  !        --> called from a spatial loop over i,j,k
  !
  subroutine calc_od_ud_vol(drho,odvol,udvol,sdetg,dx)
    real(c_double), intent(in) :: drho, sdetg, dx
    real(c_double), intent(inout) :: odvol, udvol

    if (drho>1.) then
       ! overdensity
       odvol = odvol + sdetg * dx**3
    elseif (drho<1.) then
       ! underdensity
       udvol = udvol + sdetg * dx**3
    endif

  end subroutine calc_od_ud_vol



end module spines
