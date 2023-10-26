module periodic
  !
  ! A module to implement periodic boundary conditions
  !
  implicit none

contains


  subroutine apply_periodic(j,jp1,jm1,nx,step)
    integer, intent(in) :: j, nx
    integer, intent(in), optional :: step
    integer, intent(inout) :: jp1, jm1
    integer :: stp
    
    if (present(step)) then
       stp = step
    else
       stp = 1 ! default to 1
    endif

    !
    ! Set (j+1), (j-1) depending on j value, implementing periodic BC's
    !     --> now dependent on "step" value
    !
    if (j==1) then
       jp1 = j + stp
       jm1 = nx - (stp-1)
    elseif (j==nx) then
       jp1 = stp
       jm1 = j - stp
    else
       jp1 = j + stp
       jm1 = j - stp
    endif

  end subroutine apply_periodic


  ! 
  ! Periodic boundaries when using fourth order derivatives 
  !
  subroutine apply_periodic_fourth(j,jp2,jm2,nx)
    integer, intent(in) :: j, nx
    integer, intent(inout) :: jp2, jm2
    integer :: nxm1
    nxm1 = nx-1

    if (j==1) then
       jp2 = j + 2
       jm2 = nx - 1
    elseif (j==2) then
       jp2 = j + 2
       jm2 = nx
    elseif (j==nxm1) then
       jp2 = 1
       jm2 = j - 2
    elseif (j==nx) then
       jp2 = 2
       jm2 = j - 2
    else
       jp2 = j + 2
       jm2 = j - 2
    endif

  end subroutine apply_periodic_fourth

end module periodic
