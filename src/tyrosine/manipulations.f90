module manipulations
  !
  ! A module to perform various tensor manipulations, and calculate
  !    tensor/matrix-related things:
  !
  !       - Get the (3,3) metric at a given (i,j,k) location
  !       - Invert a (3,3) matrix
  !       - Calculate the trace of a (3,3) tensor using the spatial metric
  !       - Raise and lower some indicies of rank-2 (3,3) or (4,4) tensors
  !
  use options, only:c_double
  use periodic, only:apply_periodic,apply_periodic_fourth
  implicit none


contains


  !
  ! a subroutine to return the metric (or K_ij) in (3,3) form at entire finite diff. stencil (fourth)
  !
  subroutine get_metric_at_stencil(ipos,jpos,kpos,nx,gij,gij_atstencil)
    integer, intent(in) :: ipos, jpos, kpos, nx !! positions to get metric at + grid size
    real(c_double), intent(in), dimension(6,nx,nx,nx) :: gij
    real(c_double), intent(out), dimension(13,3,3) :: gij_atstencil
    integer :: im2,im1,ip1,ip2,jm2,jm1,jp1,jp2,km2,km1,kp1,kp2

    ! -----------------------------------------------------------------------------------
    !
    ! NOTE: gij_atstencil(13,3,3) first dimension is:
    !        --> (ipos,jpos,kpos),im2,im1,ip1,ip2,jm2,jm1,jp1,jp2,km2,km1,kp1,kp2
    !
    !  and those with, e.g., im2 are (im2,j,k) i.e. other indices are at current position
    !
    ! -----------------------------------------------------------------------------------

    !
    ! apply periodic boundaries -- for whole stencil here
    !
    call apply_periodic(ipos,ip1,im1,nx)
    call apply_periodic(jpos,jp1,jm1,nx)
    call apply_periodic(kpos,kp1,km1,nx)
    call apply_periodic_fourth(ipos,ip2,im2,nx)
    call apply_periodic_fourth(jpos,jp2,jm2,nx)
    call apply_periodic_fourth(kpos,kp2,km2,nx)

    !
    ! get the metric / K_ij at all stencil positions
    !
    call get_metric_at_pos(ipos,jpos,kpos,nx,gij,gij_atstencil(1,:,:))

    call get_metric_at_pos(im2,jpos,kpos,nx,gij,gij_atstencil(2,:,:))
    call get_metric_at_pos(im1,jpos,kpos,nx,gij,gij_atstencil(3,:,:))
    call get_metric_at_pos(ip1,jpos,kpos,nx,gij,gij_atstencil(4,:,:))
    call get_metric_at_pos(ip2,jpos,kpos,nx,gij,gij_atstencil(5,:,:))

    call get_metric_at_pos(ipos,jm2,kpos,nx,gij,gij_atstencil(6,:,:))
    call get_metric_at_pos(ipos,jm1,kpos,nx,gij,gij_atstencil(7,:,:))
    call get_metric_at_pos(ipos,jp1,kpos,nx,gij,gij_atstencil(8,:,:))
    call get_metric_at_pos(ipos,jp2,kpos,nx,gij,gij_atstencil(9,:,:))

    call get_metric_at_pos(ipos,jpos,km2,nx,gij,gij_atstencil(10,:,:))
    call get_metric_at_pos(ipos,jpos,km1,nx,gij,gij_atstencil(11,:,:))
    call get_metric_at_pos(ipos,jpos,kp1,nx,gij,gij_atstencil(12,:,:))
    call get_metric_at_pos(ipos,jpos,kp2,nx,gij,gij_atstencil(13,:,:))


  end subroutine get_metric_at_stencil


  !
  ! subroutine to return the metric (3,3) at a given i,j,k position in space
  !
  subroutine get_metric_at_pos(ipos,jpos,kpos,nx,gij,gij_atpos)
    integer, intent(in) :: ipos, jpos, kpos, nx !! positions to get metric at + grid size
    real(c_double), intent(in), dimension(6,nx,nx,nx) :: gij
    real(c_double), intent(out), dimension(3,3) :: gij_atpos

    ! gij = (gxx,gxy,gxz,gyy,gyz,gzz)
    gij_atpos(1,1) = gij(1,ipos,jpos,kpos)
    gij_atpos(1,2) = gij(2,ipos,jpos,kpos)
    gij_atpos(1,3) = gij(3,ipos,jpos,kpos)
    gij_atpos(2,1) = gij(2,ipos,jpos,kpos)
    gij_atpos(2,2) = gij(4,ipos,jpos,kpos)
    gij_atpos(2,3) = gij(5,ipos,jpos,kpos)
    gij_atpos(3,1) = gij(3,ipos,jpos,kpos)
    gij_atpos(3,2) = gij(5,ipos,jpos,kpos)
    gij_atpos(3,3) = gij(6,ipos,jpos,kpos)

  end subroutine get_metric_at_pos



  !
  ! Calculate the inverse of a (3,3) matrix and output determinant
  !
  pure subroutine inv3x3(A,B,det)
    real(c_double), intent(in), dimension(3,3) :: A
    real(c_double), intent(out), dimension(3,3) :: B ! inverse matrix
    real(c_double), intent(out) :: det

    det = A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) - &
         & A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) + &
         & A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1)

    B(1,1) = A(2,2)*A(3,3) - A(2,3)*A(3,2)
    B(2,1) = A(2,3)*A(3,1) - A(2,1)*A(3,3)
    B(3,1) = A(2,1)*A(3,2) - A(2,2)*A(3,1)
    B(1,2) = A(1,3)*A(3,2) - A(1,2)*A(3,3)
    B(2,2) = A(1,1)*A(3,3) - A(1,3)*A(3,1)
    B(3,2) = A(1,2)*A(3,1) - A(1,1)*A(3,2)
    B(1,3) = A(1,2)*A(2,3) - A(1,3)*A(2,2)
    B(2,3) = A(1,3)*A(2,1) - A(1,1)*A(2,3)
    B(3,3) = A(1,1)*A(2,2) - A(1,2)*A(2,1)

    B(:,:) = B(:,:)/det
  end subroutine inv3x3


  !
  ! Calculate the trace of a (3,3) tensor given the *full* metric at all positions
  !        & position in space to calculate trace
  !
  real(c_double) function trace(A,gdown)
    !integer, intent(in) :: i,j,k,nx    ! position in space + grid size
    real(c_double), intent(in) :: A(6)    ! (3,3) symmetric tensor components
    real(c_double), intent(in) :: gdown(3,3) !,gij(6,nx,nx,nx)
    real(c_double) :: gup(3,3), detg

    !call get_metric_at_pos(i,j,k,nx,gij,gdown)
    call inv3x3(gdown,gup,detg)

    ! A = (Axx, Axy, Axz, Ayy, Ayz, Azz)
    trace = gup(1,1) * A(1) + 2.d0 * gup(1,2) * A(2) + 2.d0 * gup(1,3) * A(3) + &
         & gup(2,2) * A(4) + 2.d0 * gup(2,3) * A(5) + gup(3,3) * A(6)
  end function trace


  !
  ! Calculate the trace of a (3,3) tensor given the metric at a position in space
  !         (i.e. skips get_metric_at_pos in above)
  !
  real(c_double) function trace_atpoint(A,gdown)
    real(c_double), intent(in) :: A(6)    ! (3,3) symmetric tensor components
    real(c_double), intent(in) :: gdown(3,3)
    real(c_double) :: gup(3,3), detg

    call inv3x3(gdown,gup,detg)

    ! A = (Axx, Axy, Axz, Ayy, Ayz, Azz)
    trace_atpoint = gup(1,1) * A(1) + 2.d0 * gup(1,2) * A(2) + 2.d0 * gup(1,3) * A(3) + &
         & gup(2,2) * A(4) + 2.d0 * gup(2,3) * A(5) + gup(3,3) * A(6)

  end function trace_atpoint


  !
  ! Calculate *a single component* of the (4,4) tensor A^{mu,nu} from A_{mu,nu}
  !
  !     --> Returns A^{mu,nu} from g^{mu,alp} g^{nu,beta} A_{alp,beta}
  !
  subroutine calc_raised_comp_4D(i,j,Add,upgij,upg00,upg0i,Auu)
    integer, intent(in) :: i,j ! two down indices of tensor to raise
    real(c_double), intent(in) :: Add(4,4)
    real(c_double), intent(in) :: upgij(3,3),upg00, upg0i(3) ! components of g^{mu,nu}
    real(c_double), intent(out) :: Auu  ! component of the A_up_up tensor

    real(c_double) :: upg(4,4)
    integer :: k,l

    !
    ! fill 4D metric from 3+1 components
    upg(1,1) = upg00          ! g^{0,0}
    upg(1,2) = upg0i(1)       ! g^{0,1}
    upg(2,1) = upg0i(1)       ! g^{1,0}
    upg(1,3) = upg0i(2)       ! g^{0,2}
    upg(3,1) = upg0i(2)       ! g^{2,0}
    upg(1,4) = upg0i(3)       ! g^{0,3}
    upg(4,1) = upg0i(3)       ! g^{3,0}
    !
    ! is there a better way to do this without this loop?
    !
    do k=1,3
       do l=1,3
          upg(k+1,l+1) = upgij(k,l)
       enddo
    enddo
    !
    ! raise indices
    Auu = 0.d0
    do k=1,4
       do l=1,4
          Auu = Auu + upg(i,k) * upg(j,l) * Add(k,l)
       enddo
    enddo

  end subroutine calc_raised_comp_4D



  !
  ! Calculate a *single component* of the (3,3) tensor A^{ij} from A_{ij}
  !
  !       --> Returns A^{i,j} from \gamma^{i,k} \gamma^{j,l} A_{k,l}
  !
  subroutine calc_raised_comp(i,j,Add,upgij,Auu)
    integer, intent(in) :: i,j ! two down indices of tensor to raise
    real(c_double), intent(in) :: Add(3,3), upgij(3,3)
    real(c_double), intent(out) :: Auu  ! component of the A_up_up tensor
    integer :: k,l
    Auu = 0.d0
    do k=1,3
       do l=1,3
          Auu = Auu + upgij(i,k) * upgij(j,l) * Add(k,l)
       enddo
    enddo

  end subroutine calc_raised_comp


  !
  ! Calculate a *single component* of the (3,3) tensor A^i_j from A_{ij}
  !
  !         --> Returns A^{iu}_{id} from \gamma^{iu,k} A_{k,id}
  !
  subroutine calc_up_down_comp(iu,id,Add,upgij,Aud)
    integer, intent(in) :: iu, id   ! up and down indices of sigma to calculate
    real(c_double), intent(in), dimension(3,3) :: Add, upgij
    real(c_double), intent(out) :: Aud    ! the updown sigma returned
    integer :: k

    Aud = 0.d0
    do k=1,3
       Aud = Aud + upgij(iu,k) * Add(k,id)
    enddo

  end subroutine calc_up_down_comp

  !
  ! Calculate the *whole* (3,3) tensor A^i_j from A_{ij}
  !
  subroutine calc_up_down(Add,upgij,Aud)
    real(c_double), intent(in), dimension(3,3) :: Add, upgij
    real(c_double), intent(out) :: Aud(3,3)
    integer :: l,m,k
    Aud = 0.d0
    do l=1,3
       do m=1,3
          ! sum loop
          do k=1,3
             Aud(l,m) = Aud(l,m) + upgij(l,k) * Add(k,m)
          enddo
       enddo
    enddo
  end subroutine calc_up_down

  !
  ! Calculate the *whole* (3,3) tensor A_i^j from A_{ij}
  !
  subroutine calc_down_up(Add,upgij,Adu)
    real(c_double), intent(in), dimension(3,3) :: Add, upgij
    real(c_double), intent(out) :: Adu(3,3)
    integer :: l,m,k
    Adu = 0.d0
    do l=1,3
       do m=1,3
          ! sum loop
          do k=1,3
             Adu(l,m) = Adu(l,m) + upgij(k,m) * Add(l,k)
          enddo
       enddo
    enddo
  end subroutine calc_down_up


  !
  ! Calculate the *whole* (n,n) tensor A^a_b from A_{ab}
  !     given g_{ab}
  !
  subroutine calc_up_down_ND(n,Add,upgab,Aud)
      integer, intent(in) :: n
      real(c_double), intent(in), dimension(n,n) :: Add,upgab
      real(c_double), intent(out) :: Aud(n,n)
      integer :: l,m,k

      Aud = 0.d0
      do l=1,n
          do m=1,n
              ! sum loop
              do k=1,n
                  Aud(l,m) = Aud(l,m) + upgab(l,k) * Add(k,m)
              enddo
          enddo
      enddo
  end subroutine calc_up_down_ND



end module manipulations
