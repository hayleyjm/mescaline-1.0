module spines
  !
  ! A module containing some calculations for:
  !    - Averaging
  !    - Calculate dt(alpha), fluid R, a_FLRW
  !
  use options, only:pi,ainit,c_double,dp
  implicit none

contains

  ! ======================================================================
  !
  !                        -------- AVERAGING ---------
  !
  ! ======================================================================

  !
  ! Calculate averages over a given radius for *one sphere* only
  !
  subroutine calc_average_within_radius(xi,xj,xk,origin,radius,detg,dx,avg,posval)
    real(c_double), intent(inout) :: avg
    real(c_double), intent(in) :: radius, dx, detg
    real(c_double), intent(in) :: posval  ! quantity to be averaged at position (xi,xj,xk)
    real(c_double), intent(in) :: xi,xj,xk
    real(c_double), intent(in) :: origin(3)
    real(c_double) :: rval2

    if (radius==0.) then
       !
       ! Calculate average over whole box
       !
       avg = avg + sqrt(detg) * posval * dx**3
    else
       !
       ! Calculate where we are in radius
       !
       rval2 = (xi - origin(1))**2 + (xj - origin(2))**2 + (xk - origin(3))**2
       if (rval2<=radius**2) then
          !
          ! We are within the radius of averaging, add to the average
          !
          avg = avg + sqrt(detg) * posval * dx**3
       endif
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

    real(c_double) :: rval2, radius2, origin(3), sdetg
    integer :: l,k
    radius2 = radius*radius
    sdetg = sqrt(detg)

    if (radius==0.) then
       !
       ! Calculate avg over whole box
       !
       do k=1,navgs ! Loop over variables
          avgs(k,:) = avgs(k,:) + sdetg * dx**3 * posvals(k)
       enddo
    else
       !
       ! Calculate where we are in radius for each sphere position
       !
       do l=1,nspheres
          origin = randorigins(:,l)
          rval2 = (xi - origin(1))**2 + (xj - origin(2))**2 + (xk - origin(3))**2

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
       enddo
    endif

  end subroutine calc_average_within_radius_manyspheres



  ! ======================================================================
  !
  !           -------- MISCELLANEOUS RANDOM CALCULATIONS ---------
  !
  ! ======================================================================


  !
  ! Helps calculate the fraction of OD vs UD regions
  !        --> called from a spatial loop over i,j,k
  !
  subroutine calc_od_ud_vol(drho,odvol,udvol,sdetg,dx)
    real(c_double), intent(in) :: drho, sdetg, dx
    real(c_double), intent(inout) :: odvol, udvol

    if (drho>1._dp) then
       ! overdensity
       odvol = odvol + sdetg * dx**3
    elseif (drho<1._dp) then
       ! underdensity
       udvol = udvol + sdetg * dx**3
    endif

  end subroutine calc_od_ud_vol



end module spines
