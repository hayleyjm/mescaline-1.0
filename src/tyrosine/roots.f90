module roots
  !
  ! A module to perform the main underlying calculations ("roots") of Mescaline
  !
  !       - Ricci tensor components
  !       - Christoffel symbols
  !       - Backreaction terms and cosmological parameters
  !       - Expansion \Theta and shear scalar \sigma^2
  !
  use manipulations, only:get_metric_at_pos,inv3x3,calc_raised_comp_4D,&
    & calc_up_down,trace
  use derivatives, only:deriv1fourth,deriv1_bckwrd2nd,&
       & deriv1_bckwrd3rd,deriv1_bckwrd4th,deriv2fourth,&
       & deriv2,deriv2_mix,deriv2_mixfourth,deriv2_bckwrd2nd,deriv1_forward4th
  use periodic, only:apply_periodic,apply_periodic_fourth
  use options, only:dp,pi,c_double,dit,tinitial,dtfac,&
    & tderivs,nord_dt,writeQdterms,domain_type
  use prints, only:print_info,print_error
  use tripreports, only:write_avg
  implicit none

contains


    !
    ! Subroutine to get Christoffel symbols and trK on the whole grid, at ONE TIME
    !     called by the raytracer if first=False
    !
    subroutine get_christoffels_trK_allgrid(nx,dx,gij,kij,gamijk,tracek)
        integer, intent(in) :: nx
        real(c_double), intent(in) :: dx
        real(c_double), intent(in), dimension(6,nx,nx,nx) :: gij,kij

        real(c_double), intent(out) :: gamijk(3,3,3,nx,nx,nx),tracek(nx,nx,nx)

        real(c_double) :: gdown(3,3),gup(3,3),detg
        integer :: i,j,k,l,m,n

        !$omp parallel do default (none) &
        !$omp shared(nx,dx,gij,kij,gamijk,tracek) &
        !$omp private(i,j,k,l,m,n,gdown)
        do k=1,nx
           do j=1,nx

              do i=1,nx
                  !
                  ! calculate tracek
                  call get_metric_at_pos(i,j,k,nx,gij,gdown)
                  tracek(i,j,k) = trace(kij(:,i,j,k),gdown)
                  !
                  ! invert gdown to pass into christofel and save number of inv3x3 calls
                  call inv3x3(gdown,gup,detg)

                  !
                  ! Loop over components of Christoffel symbols
                  do l=1,3
                      do m=1,3
                          do n=1,3
                              call get_christoffel(i,j,k,gij,nx,dx,l,m,n,&
                              & gamijk(l,m,n,i,j,k),gup)
                          enddo
                      enddo
                  enddo

                enddo
            enddo
        enddo
        !$omp end parallel do

    end subroutine get_christoffels_trK_allgrid



    !
    ! Subroutine to get trK on the whole grid, at ONE TIME
    !     called by the raytracer (testing speed w/o Christoffels)
    !
    subroutine get_trK_allgrid(nx,dx,gij,kij,tracek)
        integer, intent(in) :: nx
        real(c_double), intent(in) :: dx
        real(c_double), intent(in), dimension(6,nx,nx,nx) :: gij,kij

        real(c_double), intent(out) :: tracek(nx,nx,nx)

        real(c_double) :: gdown(3,3)
        integer :: i,j,k

        !$omp parallel do default (none) &
        !$omp shared(nx,dx,gij,kij,tracek) &
        !$omp private(i,j,k,gdown)
        do k=1,nx
           do j=1,nx

              do i=1,nx
                  !
                  ! calculate tracek
                  call get_metric_at_pos(i,j,k,nx,gij,gdown)
                  tracek(i,j,k) = trace(kij(:,i,j,k),gdown)

                enddo
            enddo
        enddo
        !$omp end parallel do

    end subroutine get_trK_allgrid




    !
    ! Subroutine to get Christoffel symbols and trK on the whole grid, at all times
    !      because we need to do this at the beginning of mescaline
    !   ** we get Gammas and K diffferently for RT now; so this is only called from mesc routine
    !     when we always have Gammas at one time only, i.e.nt_gam=1.
    !
    subroutine get_christoffels_trK_allgrid_alltime(nx,nt,ti,dx,gij,kij,gamijk,tracek)
        integer, intent(in) :: nx,nt,ti
        real(c_double), intent(in) :: dx
        real(c_double), intent(in), dimension(6,nx,nx,nx,nt) :: gij,kij

        real(c_double), intent(out) :: gamijk(3,3,3,nx,nx,nx),tracek(nx,nx,nx,nt)

        real(c_double) :: gdown(3,3)
        integer :: i,j,k,l,m,n

        !$omp parallel do default (none) &
        !$omp shared(nx,nt,ti,dx,gij,kij,gamijk,tracek) &
        !$omp private(i,j,k,l,m,n,gdown)
        do k=1,nx
           do j=1,nx

              do i=1,nx
                  !
                  ! Loop over components of Christoffel symbols
                  do l=1,3
                      do m=1,3
                          do n=1,3
                              ! Only get these at one time; this routine only called from ricci
                              call get_christoffel(i,j,k,gij(:,:,:,:,ti),nx,dx,l,m,n,&
                              & gamijk(l,m,n,i,j,k))
                          enddo
                      enddo
                  enddo
                  !
                  ! calculate tracek: loop over time and get whole time stencil
                  do n=1,nt
                      call get_metric_at_pos(i,j,k,nx,gij(:,:,:,:,n),gdown)
                      tracek(i,j,k,n) = trace(kij(:,i,j,k,n),gdown)
                  enddo

                enddo
            enddo
        enddo
        !$omp end parallel do

    end subroutine get_christoffels_trK_allgrid_alltime



  !
  ! subroutine to calculate a given component of the (3,3) Ricci tensor
  !
  subroutine get_ricci_component(i,j,k,nx,dx,gij,ridx1,ridx2,riccic,Chrsijk)
    integer, intent(in) :: i,j,k,nx             ! current position in space & grid size
    real(c_double), intent(in) :: dx
    ! Christoffel full grid array is optional; if not present, calculate them
    !    --> this is for the raytracer, because it's not necessarily better to calculate all of these beforehand
    real(c_double), intent(in) :: gij(6,nx,nx,nx)
    real(c_double), intent(in), optional :: Chrsijk(3,3,3,nx,nx,nx)
    integer, intent(in) :: ridx1, ridx2         ! the indices of the component, R_{idx1,idx2}

    real(c_double), intent(out) :: riccic
    real(c_double) :: term1, term2, term3, term4

    integer :: l, m                             ! dummy indices of \Gamma's
    real(c_double) :: gamma1_p1, gamma1_m1, gamma2_p1, gamma2_m1, gamma3a, gamma3b, gamma4a, gamma4b
    real(c_double) :: gamma1_p2, gamma1_m2, gamma2_p2, gamma2_m2
    integer :: ip1,im1,jp1,jm1,kp1,km1,ip2,im2,jp2,jm2,kp2,km2 ! for periodic boundaries
    !
    ! A note on notation: here gamma1 indicates the christoffel symbol in the first term of the Ricci tensor (etc...):
    ! R_{idx1,idx2} = \partial_{i}\Gamma^{i}_{idx1,idx2} - \partial_{idx2}\Gamma^{l}_{idx1,l}
    ! + \Gamma^{l}_{n,l}\Gamma^{n}_{idx1,idx2} - \Gamma^{l}_{idx2,n}\Gamma^{n}_{idx1,l}
    !

    ! Apply periodic boundary conditions
    call apply_periodic(i,ip1,im1,nx)
    call apply_periodic(j,jp1,jm1,nx)
    call apply_periodic(k,kp1,km1,nx)
    call apply_periodic_fourth(i,ip2,im2,nx)
    call apply_periodic_fourth(j,jp2,jm2,nx)
    call apply_periodic_fourth(k,kp2,km2,nx)
    !
    ! initialise sums to zero
    term1 = 0._dp
    term2 = 0._dp
    term3 = 0._dp
    term4 = 0._dp

    if (present(Chrsijk)) then
        ! we have Christoffels in an array; use them
        do l=1,3
            select case (l)
            case(1)
                ! x-deriv: calculate \Gamma1 at i+1, i-1
                !
                ! FIRST TERM
                gamma1_p1 = Chrsijk(l,ridx1,ridx2,ip1,j,k)
                gamma1_m1 = Chrsijk(l,ridx1,ridx2,im1,j,k)
                gamma1_p2 = Chrsijk(l,ridx1,ridx2,ip2,j,k)
                gamma1_m2 = Chrsijk(l,ridx1,ridx2,im2,j,k)

            case(2)
                ! y-deriv: calculate \Gamma1 at j+1, j-1
                !
                ! FIRST TERM
                gamma1_p1 = Chrsijk(l,ridx1,ridx2,i,jp1,k)
                gamma1_m1 = Chrsijk(l,ridx1,ridx2,i,jm1,k)
                gamma1_p2 = Chrsijk(l,ridx1,ridx2,i,jp2,k)
                gamma1_m2 = Chrsijk(l,ridx1,ridx2,i,jm2,k)

            case(3)
                ! z-deriv: calculate \Gamma1 at k+1, k-1
                !
                ! FIRST TERM
                gamma1_p1 = Chrsijk(l,ridx1,ridx2,i,j,kp1)
                gamma1_m1 = Chrsijk(l,ridx1,ridx2,i,j,km1)
                gamma1_p2 = Chrsijk(l,ridx1,ridx2,i,j,kp2)
                gamma1_m2 = Chrsijk(l,ridx1,ridx2,i,j,km2)
            end select
            term1 = term1 + deriv1fourth(gamma1_p1,gamma1_p2,gamma1_m1,gamma1_m2,dx)

           !
           ! SECOND TERM
           select case (ridx2)
           case(1)
               !
               ! take the x-derivative of christoffel, get at i-stencil
               gamma2_p1 = Chrsijk(l,l,ridx1,ip1,j,k)
               gamma2_m1 = Chrsijk(l,l,ridx1,im1,j,k)
               gamma2_p2 = Chrsijk(l,l,ridx1,ip2,j,k)
               gamma2_m2 = Chrsijk(l,l,ridx1,im2,j,k)

           case(2)
               !
               ! take the y-derivaitve of christoffel, get at j-stencil
               gamma2_p1 = Chrsijk(l,l,ridx1,i,jp1,k)
               gamma2_m1 = Chrsijk(l,l,ridx1,i,jm1,k)
               gamma2_p2 = Chrsijk(l,l,ridx1,i,jp2,k)
               gamma2_m2 = Chrsijk(l,l,ridx1,i,jm2,k)

           case(3)
               !
               ! take the z-derivaitve of christoffel, get at j-stencil
               gamma2_p1 = Chrsijk(l,l,ridx1,i,j,kp1)
               gamma2_m1 = Chrsijk(l,l,ridx1,i,j,km1)
               gamma2_p2 = Chrsijk(l,l,ridx1,i,j,kp2)
               gamma2_m2 = Chrsijk(l,l,ridx1,i,j,km2)

            end select
            term2 = term2 + deriv1fourth(gamma2_p1,gamma2_p2,gamma2_m1,gamma2_m2,dx)

           ! additional sum for last two terms (\Gamma*\Gamma terms)
           do m=1,3
              !
              ! THIRD TERM
              gamma3a = Chrsijk(l,m,l,i,j,k)
              gamma3b = Chrsijk(m,ridx1,ridx2,i,j,k)

              term3 = term3 + gamma3a * gamma3b

              !
              ! FOURTH TERM
              gamma4a = Chrsijk(l,ridx2,m,i,j,k)
              gamma4b = Chrsijk(m,ridx1,l,i,j,k)

              term4 = term4 + gamma4a * gamma4b
           enddo
        enddo

    else
        ! We don't have Christoffels as an array; need to calculate them!
        do l=1,3
            select case (l)
            case(1)
                ! x-deriv: calculate \Gamma1 at i+1, i-1
                !
                ! FIRST TERM
                call get_christoffel(ip1,j,k,gij,nx,dx,l,ridx1,ridx2,gamma1_p1)
                call get_christoffel(im1,j,k,gij,nx,dx,l,ridx1,ridx2,gamma1_m1)
                call get_christoffel(ip2,j,k,gij,nx,dx,l,ridx1,ridx2,gamma1_p2)
                call get_christoffel(im2,j,k,gij,nx,dx,l,ridx1,ridx2,gamma1_m2)

            case(2)
                ! y-deriv: calculate \Gamma1 at j+1, j-1
                !
                ! FIRST TERM
                call get_christoffel(i,jp1,k,gij,nx,dx,l,ridx1,ridx2,gamma1_p1)
                call get_christoffel(i,jm1,k,gij,nx,dx,l,ridx1,ridx2,gamma1_m1)
                call get_christoffel(i,jp2,k,gij,nx,dx,l,ridx1,ridx2,gamma1_p2)
                call get_christoffel(i,jm2,k,gij,nx,dx,l,ridx1,ridx2,gamma1_m2)

            case(3)
                ! z-deriv: calculate \Gamma1 at k+1, k-1
                !
                ! FIRST TERM
                call get_christoffel(i,j,kp1,gij,nx,dx,l,ridx1,ridx2,gamma1_p1)
                call get_christoffel(i,j,km1,gij,nx,dx,l,ridx1,ridx2,gamma1_m1)
                call get_christoffel(i,j,kp2,gij,nx,dx,l,ridx1,ridx2,gamma1_p2)
                call get_christoffel(i,j,km2,gij,nx,dx,l,ridx1,ridx2,gamma1_m2)

            end select
            term1 = term1 + deriv1fourth(gamma1_p1,gamma1_p2,gamma1_m1,gamma1_m2,dx)

           !
           ! SECOND TERM
           select case (ridx2)
           case(1)
               !
               ! take the x-derivative of christoffel, get at i-stencil
               call get_christoffel(ip1,j,k,gij,nx,dx,l,l,ridx1,gamma2_p1)
               call get_christoffel(im1,j,k,gij,nx,dx,l,l,ridx1,gamma2_m1)
               call get_christoffel(ip2,j,k,gij,nx,dx,l,l,ridx1,gamma2_p2)
               call get_christoffel(im2,j,k,gij,nx,dx,l,l,ridx1,gamma2_m2)

           case(2)
               !
               ! take the y-derivaitve of christoffel, get at j-stencil
               call get_christoffel(i,jp1,k,gij,nx,dx,l,l,ridx1,gamma2_p1)
               call get_christoffel(i,jm1,k,gij,nx,dx,l,l,ridx1,gamma2_m1)
               call get_christoffel(i,jp2,k,gij,nx,dx,l,l,ridx1,gamma2_p2)
               call get_christoffel(i,jm2,k,gij,nx,dx,l,l,ridx1,gamma2_m2)

           case(3)
               !
               ! take the z-derivaitve of christoffel, get at j-stencil
               call get_christoffel(i,j,kp1,gij,nx,dx,l,l,ridx1,gamma2_p1)
               call get_christoffel(i,j,km1,gij,nx,dx,l,l,ridx1,gamma2_m1)
               call get_christoffel(i,j,kp2,gij,nx,dx,l,l,ridx1,gamma2_p2)
               call get_christoffel(i,j,km2,gij,nx,dx,l,l,ridx1,gamma2_m2)

            end select
            term2 = term2 + deriv1fourth(gamma2_p1,gamma2_p2,gamma2_m1,gamma2_m2,dx)

           ! additional sum for last two terms (\Gamma*\Gamma terms)
           do m=1,3
              !
              ! THIRD TERM
              call get_christoffel(i,j,k,gij,nx,dx,l,m,l,gamma3a)
              call get_christoffel(i,j,k,gij,nx,dx,m,ridx1,ridx2,gamma3b)

              term3 = term3 + gamma3a * gamma3b
              !
              ! FOURTH TERM
              call get_christoffel(i,j,k,gij,nx,dx,l,ridx2,m,gamma4a)
              call get_christoffel(i,j,k,gij,nx,dx,m,ridx1,l,gamma4b)
              term4 = term4 + gamma4a * gamma4b
           enddo
        enddo

    endif

    riccic = term1 - term2 + term3 - term4
  end subroutine get_ricci_component






  !
  ! a subroutine to take in current poisition in space, and indices of the christoffel and calc christoffel
  !
  !  --> * checked March 27th, 2020 (Hi from lockdown!) and all OK
  !
  subroutine get_christoffel(i,j,k,gij,nx,dx,idx1,idx2,idx3,chr,gupijk)
    integer, intent(in) :: nx                        ! the size of the grid (assuming uniform)
    integer, intent(in) :: i,j,k                     ! current position in space
    integer, intent(in) :: idx1, idx2, idx3          ! the indices of the christoffel: \Gamma^{idx1}_{idx2,idx3}
    real(c_double), intent(in), dimension(6,nx,nx,nx) :: gij
    real(c_double), intent(in) :: dx                       ! the gridspacing
    real(c_double), intent(in), optional :: gupijk(3,3) ! g^ij at i,j,k -- optional because sometimes we have it, sometimes we don't
    real(c_double), intent(out) :: chr                     ! the resulting christoffel symbol
    real(c_double), dimension(3) :: g3dum_p1, g3dum_m1, gdum2_m1, gdum2_p1, upg1dum  ! the g_{k,l} and g_{l,j} terms we pass in depending on which derivative we are taking of these terms, see explanation of notation below...
    real(c_double), dimension(3) :: g3dum_p2, g3dum_m2, gdum2_p2, gdum2_m2
    real(c_double) :: g23_ip1, g23_im1, g23_jp1, g23_jm1, g23_kp1, g23_km1  ! more temp gij terms
    real(c_double) :: g23_ip2, g23_im2, g23_jp2, g23_jm2, g23_kp2, g23_km2  ! more temp gij terms
    real(c_double) :: g23p1(3), g23m1(3), detg, g23p2(3), g23m2(3)
    real(c_double), dimension(3,3) :: gup, gdown, gdown_ip1, gdown_im1, gdown_jp1, gdown_jm1, gdown_kp1, gdown_km1
    real(c_double), dimension(3,3) :: gdown_ip2,gdown_im2,gdown_jp2,gdown_jm2,gdown_kp2,gdown_km2
    real(c_double) :: trm1, trm2, trm3
    integer :: l, ip1,im1,jp1,jm1,kp1,km1,ip2,im2,jp2,jm2,kp2,km2

    ! Apply periodic boundary conditions
    call apply_periodic(i,ip1,im1,nx)
    call apply_periodic(j,jp1,jm1,nx)
    call apply_periodic(k,kp1,km1,nx)
    call apply_periodic_fourth(i,ip2,im2,nx)
    call apply_periodic_fourth(j,jp2,jm2,nx)
    call apply_periodic_fourth(k,kp2,km2,nx)

    ! fill 3-metric's with values at certain positions
    !
    call get_metric_at_pos(i,j,k,nx,gij,gdown)
    call get_metric_at_pos(ip1,j,k,nx,gij,gdown_ip1)
    call get_metric_at_pos(im1,j,k,nx,gij,gdown_im1)
    call get_metric_at_pos(i,jp1,k,nx,gij,gdown_jp1)
    call get_metric_at_pos(i,jm1,k,nx,gij,gdown_jm1)
    call get_metric_at_pos(i,j,kp1,nx,gij,gdown_kp1)
    call get_metric_at_pos(i,j,km1,nx,gij,gdown_km1)
    call get_metric_at_pos(ip2,j,k,nx,gij,gdown_ip2)
    call get_metric_at_pos(im2,j,k,nx,gij,gdown_im2)
    call get_metric_at_pos(i,jp2,k,nx,gij,gdown_jp2)
    call get_metric_at_pos(i,jm2,k,nx,gij,gdown_jm2)
    call get_metric_at_pos(i,j,kp2,nx,gij,gdown_kp2)
    call get_metric_at_pos(i,j,km2,nx,gij,gdown_km2)
    if (present(gupijk)) then
        ! we've been passed in g^ij -- don't invert!
        gup = gupijk
    else
        ! we haven't been passed in g^ij -- invert!
        call inv3x3(gdown,gup,detg)
    endif
    !
    ! Notes on notation below: g3dum is: g_{idx3,dummy}, gdum2 is: g_{dummy,idx2}
    ! where idx2, idx3 are the indices of the christoffel: \Gamma^{idx1}_{idx2, idx3}
    ! upg1dum is: g^{idx1,dummy}, g23 is: g_{idx2, idx3}
    !

    ! extract parts of the metric to send in to be calculated
    select case(idx2)
    case(1)
       ! we are taking the x-deriv of g_{k,l}
       g3dum_p1 = gdown_ip1(idx3,:)
       g3dum_m1 = gdown_im1(idx3,:)
       g3dum_p2 = gdown_ip2(idx3,:)
       g3dum_m2 = gdown_im2(idx3,:)
    case(2)
       ! we are taking the y-deriv of g_{k,l}
       g3dum_p1 = gdown_jp1(idx3,:)
       g3dum_m1 = gdown_jm1(idx3,:)
       g3dum_p2 = gdown_jp2(idx3,:)
       g3dum_m2 = gdown_jm2(idx3,:)
    case(3)
       ! we are taking the z-deriv of g_{k,l}
       g3dum_p1 = gdown_kp1(idx3,:)
       g3dum_m1 = gdown_km1(idx3,:)
       g3dum_p2 = gdown_kp2(idx3,:)
       g3dum_m2 = gdown_km2(idx3,:)
    end select

    select case(idx3)
    case(1)
       ! we are taking the x-deriv of g_{l,j}
       gdum2_p1 = gdown_ip1(:,idx2)
       gdum2_m1 = gdown_im1(:,idx2)
       gdum2_p2 = gdown_ip2(:,idx2)
       gdum2_m2 = gdown_im2(:,idx2)
     case(2)
       ! we are taking the y-deriv of g_{l,j}
       gdum2_p1 = gdown_jp1(:,idx2)
       gdum2_m1 = gdown_jm1(:,idx2)
       gdum2_p2 = gdown_jp2(:,idx2)
       gdum2_m2 = gdown_jm2(:,idx2)
    case(3)
       ! we are taking the z-deriv of g_{l,j}
       gdum2_p1 = gdown_kp1(:,idx2)
       gdum2_m1 = gdown_km1(:,idx2)
       gdum2_p2 = gdown_kp2(:,idx2)
       gdum2_m2 = gdown_km2(:,idx2)
    end select

    upg1dum = gup(idx1,:)
    g23_ip1 = gdown_ip1(idx2,idx3)
    g23_im1 = gdown_im1(idx2,idx3)
    g23_jp1 = gdown_jp1(idx2,idx3)
    g23_jm1 = gdown_jm1(idx2,idx3)
    g23_kp1 = gdown_kp1(idx2,idx3)
    g23_km1 = gdown_km1(idx2,idx3)
    g23_ip2 = gdown_ip2(idx2,idx3)
    g23_im2 = gdown_im2(idx2,idx3)
    g23_jp2 = gdown_jp2(idx2,idx3)
    g23_jm2 = gdown_jm2(idx2,idx3)
    g23_kp2 = gdown_kp2(idx2,idx3)
    g23_km2 = gdown_km2(idx2,idx3)

    ! for the term g_{idx2,idx3} put into list to take certain derivs in loop
    g23p1 = (/ g23_ip1, g23_jp1, g23_kp1 /)
    g23m1 = (/ g23_im1, g23_jm1, g23_km1 /)
    g23p2 = (/ g23_ip2, g23_jp2, g23_kp2 /)
    g23m2 = (/ g23_im2, g23_jm2, g23_km2 /)

    chr = 0._dp
    do l=1,3
        trm1 = deriv1fourth(g3dum_p1(l),g3dum_p2(l),g3dum_m1(l),g3dum_m2(l),dx)
        trm2 = deriv1fourth(gdum2_p1(l),gdum2_p2(l),gdum2_m1(l),gdum2_m2(l),dx)
        trm3 = deriv1fourth(g23p1(l),g23p2(l),g23m1(l),g23m2(l),dx)
       chr = chr + 0.5_dp * upg1dum(l) * (trm1 + trm2 - trm3)  ! summation
    enddo

  end subroutine get_christoffel




  !
  !    subroutine to calculate backreaction terms and omegas from avgs
  !  -- uses Buchert+ 2018's new formalism via private communication --
  !
  subroutine get_backreaction_omegas(it,time,rad,nspheres,randorigins,avgs,lenavgs,rhoavgall,rhoavgall_notilde,&
       & vDball,Qd,hub,omegam,omegaR,omegaQ,delta,aDb)
    integer, intent(in) :: nspheres,lenavgs,it ! the number of average values we have
    real(c_double), intent(in)    :: time,rad,vDball
    real(c_double), intent(in), dimension(3,nspheres) :: randorigins
    real(c_double), intent(inout) :: avgs(lenavgs,nspheres),rhoavgall,rhoavgall_notilde
    real(c_double), intent(out), dimension(nspheres) :: Qd,hub,omegam,omegaR,omegaQ,delta
    real(c_double), intent(out), dimension(nspheres) :: aDb ! effective scale factor
    character(len=100) :: message,loc
    loc = " get_backreaction_omegas" ! our routine location

    !
    !    -- normalise averages by sum of lorentz factor --
    !      this simplifies from the relation between the
    !      averages w.r.t b and h. see "user guide"
    !
    avgs(1,:) = avgs(1,:) / avgs(7,:) ! \int_D( \Gamma \tilde{\Theta} ) / \int_D \Gamma
    avgs(2,:) = avgs(2,:) / avgs(7,:) ! \int_D( \Gamma \tilde{R} )      / \int_D \Gamma
    avgs(3,:) = avgs(3,:) / avgs(7,:) ! \int_D( \Gamma \tilde{\rho} )   / \int_D \Gamma
    avgs(4,:) = avgs(4,:) / avgs(7,:) ! \int_D( \Gamma \tilde{\sigma^2} )    / \int_D \Gamma
    avgs(5,:) = avgs(5,:) / avgs(7,:) ! \int_D( \Gamma \tilde{\Theta^2} )    / \int_D \Gamma
    avgs(6,:) = avgs(6,:) / avgs(7,:) ! \int_D( \Gamma \tilde{vorticity^2} ) / \int_D \Gamma

    !
    ! calculate kinematical backreaction - see Buchert+2018 private communication
    !
    Qd = (2._dp / 3._dp) * ( avgs(5,:) - avgs(1,:)**2 ) - 2._dp * avgs(4,:) + 2._dp * avgs(6,:)

    ! write Qd terms if we want them
    if (writeQdterms) then
        call write_avg(avgs(5,:),"theta2_avg",it,time,nspheres,rad,domain_type,&
            & "\eta, <\Theta^2>",randorigins)
        call write_avg(avgs(4,:),"sigma2_avg",it,time,nspheres,rad,domain_type,&
            & "\eta, <\sigma^2>",randorigins)
        call write_avg(avgs(6,:),"omega2_avg",it,time,nspheres,rad,domain_type,&
            & "\eta, <\omega^2>",randorigins)
    endif

    !
    ! calculate omegas
    !
    hub = avgs(1,:) / 3._dp
    omegaR = - avgs(2,:) / (6._dp * hub*hub)
    omegam = 8._dp * pi * avgs(3,:) / (3._dp * hub*hub)
    omegaQ = - Qd / (6._dp * hub*hub)

    if (nspheres==1) then
       write(message,"(a,ES15.8)") "          Qd = ",Qd
       call print_info(message,loc)
       write(message,"(a,ES15.8)") "      omegaQ = ",omegaQ
       call print_info(message,loc)
       write(message,"(a,ES15.8)") "      omegaR = ",omegaR
       call print_info(message,loc)
       write(message,"(a,ES15.8)") "      omegam = ",omegam
       call print_info(message,loc)
       write(message,"(a,ES15.8)") "         hub = ",hub
       call print_info(message,loc)
    else
       write(message,"(a,ES15.8)") "     mean Qd = ",sum(Qd)/nspheres
       call print_info(message,loc)
       write(message,"(a,ES15.8)") " mean omegaQ = ",sum(omegaQ)/nspheres
       call print_info(message,loc)
       write(message,"(a,ES15.8)") " mean omegaR = ",sum(omegaR)/nspheres
       call print_info(message,loc)
       write(message,"(a,ES15.8)") " mean omegam = ",sum(omegam)/nspheres
       call print_info(message,loc)
       write(message,"(a,ES15.8)") "    mean hub = ",sum(hub)/nspheres
       call print_info(message,loc)
    endif

    !
    ! calculate delta for each sphere
    rhoavgall         = rhoavgall / vDball         ! global average density
    rhoavgall_notilde = rhoavgall_notilde / vDball ! global average density -- NO tilde scaling (for dL calc)
    delta     = avgs(3,:) / rhoavgall - 1._dp      ! note avgs(3,:) HAS tilde scaling, so we want it here too
    !
    ! calculate effective scale factor expansion using V^b_D
    !
    ! aDb = (V^b_D / V^b_D(t0))^(1/3) = (V_D * <Gamma> / V^b_D(t0))^(1/3)
    !                                 = ( \int_D \Gamma / V^b_D(t0) )^(1/3)
    !
    ! aDb = ( avgs(7,:) / vDbt0 )**(1._dp / 3._dp)
    !
    ! We remove the vDbt0 scaling so we are consistent when ainit/=1
    !      this still gives aDb~1 (incl. perturbs) if you have ainit=1
    aDb = ( avgs(7,:) )**(1._dp / 3._dp)
    if (rad==0.) then
       write(message,"(a,ES12.5)") "          aD = ",aDb
       call print_info(message,loc)
    endif

  end subroutine get_backreaction_omegas


  !
  ! A subroutine to return the 3-d fluid-restframe curvature as defined in Buchert+2019
  !
  ! curlyR := \nabla_mu u^nu \nabla_nu u^mu - \Theta^2 + R + 2R_munu u^mu u^nu
  !
  !  (here, R is 4Ricci scalar, R_munu is 4 Ricci tensor)
  !
  subroutine get_fluid_curvature(nx,nt,dx,dt,ipos,jpos,kpos,gdown_atstencil,gup,kij,&
       & kdown_atstencil,tracek,alp,threeRij,gamijk,theta,lorentz,fvelU_atstencil,dtudt,dtudi,&
       & dtalp,dialp,diudt,kud_atstencil,aUmu,admu,maga,fluidR,fourR,fourRdd,Rvsq)
    integer, intent(in) :: nx,nt,ipos,jpos,kpos
    real(c_double), intent(in) :: dx,dt
    real(c_double), intent(in) :: theta,lorentz,fvelU_atstencil(13,3) ! u^i
    real(c_double), intent(in) :: dialp(3),dtalp,dtudt,dtudi(3),diudt(3)
    real(c_double), intent(in) :: gup(3,3),threeRij(3,3),gamijk(3,3,3) ! 3R_ij and Christoffels
    real(c_double), intent(in) :: kij(6,nx,nx,nx,nt),tracek(nx,nx,nx,nt)
    real(c_double), dimension(13,3,3),   intent(in) :: gdown_atstencil,kdown_atstencil
    real(c_double), dimension(nx,nx,nx), intent(in) :: alp

    ! output 3D fluid R
    real(c_double), intent(out) :: fluidR,fourR,fourRdd(4,4),Rvsq
    real(c_double), intent(out) :: aUmu(4),admu(4),maga         ! four-acceleration components and magnitude
    real(c_double), dimension(13,3,3), intent(out) :: kud_atstencil ! for violation_inloop

    real(c_double), dimension(3,3) :: fourRij,kdown,kud,diuUj,gdown
    real(c_double), dimension(3)   :: aUi,adi,fvelU,fvelUip2,fvelUip1,fvelUim1,fvelUim2,fourR0i
    real(c_double), dimension(3)   :: fvelUjp2,fvelUjp1,fvelUjm1,fvelUjm2,fvelUkp2,fvelUkp1,fvelUkm1,fvelUkm2
    real(c_double) :: dtuUi(3),diuUt(3),uUt,aU0!,adi(3)
    real(c_double) :: alpijk,RUUsum,nabu_nabu,fourR00!,fourR,
    real(c_double) :: term1,term2a(3),term2b(3),term3a(3,3),term3b(3,3)
    !real(c_double) :: at,rhot,asq,rhostar,flrwhub,adot,addot

    integer :: i,j,k
    !
    ! DOUBLE 'nt' in call below is 'nt,ti' -- ti is index in time dimension of 'current' time
    !     (since 4D ricci routine called from RT where this is ti=1)
    call get_4D_ricci_tensor(nx,nt,nt,dx,dt,ipos,jpos,kpos,gdown_atstencil,gup,kij,&
         & kdown_atstencil,tracek,alp,threeRij,gamijk,dialp,kud_atstencil,fourR00,&
         & fourR0i,fourRij,fourR)
    ! Store xx and xy components for test
    fourRdd(1,1)   = fourR00
    fourRdd(1,2:4) = fourR0i
    fourRdd(2:4,1) = fourR0i
    fourRdd(2:4,2:4) = fourRij
    !print*, 'R00, R0i = ',fourR00,fourR0i
    !
    ! store K^i_j, K_ij (and some other things to limit memory access)
    kud    = kud_atstencil(1,:,:)
    kdown  = kdown_atstencil(1,:,:)
    gdown  = gdown_atstencil(1,:,:)
    alpijk = alp(ipos,jpos,kpos)
    uUt    = lorentz / alpijk      ! u^0 four-vel

    ! ----------------------------------------
    !
    ! calc \pd_i u^j
    !
    ! *_atstencil first dim is:
    ! --> (pos,im2,im1,ip1,ip2,jm2,jm1,jp1,jp2,km2,km1,kp1,kp2)
    !
    fvelU = fvelU_atstencil(1,:)
    fvelUim2 = fvelU_atstencil(2,:)
    fvelUim1 = fvelU_atstencil(3,:)
    fvelUip1 = fvelU_atstencil(4,:)
    fvelUip2 = fvelU_atstencil(5,:)
    fvelUjm2 = fvelU_atstencil(6,:)
    fvelUjm1 = fvelU_atstencil(7,:)
    fvelUjp1 = fvelU_atstencil(8,:)
    fvelUjp2 = fvelU_atstencil(9,:)
    fvelUkm2 = fvelU_atstencil(10,:)
    fvelUkm1 = fvelU_atstencil(11,:)
    fvelUkp1 = fvelU_atstencil(12,:)
    fvelUkp2 = fvelU_atstencil(13,:)

    !
    ! initialise R_ij v^i v^j = R_ij u^i u^j / \Gam^2 sum
    Rvsq = 0._dp
    do j=1,3
       do i=1,3

          Rvsq = Rvsq + threeRij(i,j) * fvelU(i) * fvelU(j) / lorentz**2

          select case(i)
          case(1)
             ! x-derivative
             diuUj(i,j) = deriv1fourth(fvelUip1(j),fvelUip2(j),fvelUim1(j),fvelUim2(j),dx)
          case(2)
             ! y-derivative
             diuUj(i,j) = deriv1fourth(fvelUjp1(j),fvelUjp2(j),fvelUjm1(j),fvelUjm2(j),dx)
          case(3)
             ! z-derivative
             diuUj(i,j) = deriv1fourth(fvelUkp1(j),fvelUkp2(j),fvelUkm1(j),fvelUkm2(j),dx)
          end select
       enddo
    enddo
    ! ----------------------------------------

    ! Initialise some sums...
    nabu_nabu = 0._dp; dtUui = 0._dp

    term1  = - (dtudt + lorentz * dtalp ) / alpijk**2
    term2a = 0._dp; term2b = 0._dp
    term3a = 0._dp; term3b = 0._dp
    diuUt  = (2._dp * lorentz * dialp - diudt) / alpijk**2
    ! ---------------------------
    !
    ! Loop breakdowns
    !    k:
    !  -- sum loop for term1
    !  -- only component loop for dtuUi
    !  -- only component loop for term2's
    !  -- 1st component loop for term3s
    !
    !  -- component loop for accel^k
    !  -- sum loop for accel^0
    !
    ! ---------------------------
    aU0 = 0._dp; aUi = 0._dp; adi = 0._dp; maga = 0._dp
    !
    ! zero-sum terms for a^0
    aU0 = - uUt * dtudt / alpijk**2 - uUt**2 * dtalp / alpijk

    do k=1,3

       term1 = term1 + fvelU(k) * dialp(k) / alpijk

       term2b(k) = diuUt(k) + lorentz * dialp(k) / alpijk**2

       ! single-sum terms for a^0  (no zero-sum for a^k)
       aU0 = aU0 + fvelU(k) * diuUt(k) + 2._dp * uUt * fvelU(k) * dialp(k) / alpijk

       ! ---------------------------
       !    j:
       !  -- sum loop for term2's
       !  -- 2nd component loop for term3's
       !
       !  -- sum loop for accel^k
       !  -- sum loop for accel^0
       ! ---------------------------
       do j=1,3

          dtuUi(k)  = dtuUi(k) + gup(k,j) * dtudi(j) + 2._dp * alpijk * fvelU(j) * kud(k,j)
          term2a(k) = term2a(k) + lorentz * gup(k,j) * dialp(j) - alpijk * kud(k,j) * fvelU(j)! + dtuUi
          term2b(k) = term2b(k) - fvelU(j) * kdown(k,j) / alpijk

          ! note indices are SUPPOSED to be swapped around here - this term is confusing. See your notes. Pretty sure this is ok.
          term3a(j,k) = diuUj(j,k) - lorentz * kud(k,j)
          term3b(j,k) = diuUj(k,j) - lorentz * kud(j,k)

          ! double-sum terms for a^0 (single-sum for a^k)
          aU0    = aU0 - kdown(j,k) * fvelU(j) * fvelU(k) / alpijk
          aUi(k) = aUi(k) + uUt * gup(j,k) * dtudi(j) + fvelU(j) * diuUj(j,k) + uUt**2 * alpijk * gup(j,k) * dialp(j)


          ! ---------------------------
          !  i:
          !  -- sum loop for term3's
          !
          !  -- sum loop for accel^k
          ! ---------------------------
          do i=1,3
             term3a(j,k) = term3a(j,k) + gamijk(k,j,i) * fvelU(i)
             term3b(j,k) = term3b(j,k) + gamijk(j,k,i) * fvelU(i)

             ! double-sum terms for a^k
             aUi(k) = aUi(k) + fvelU(i) * fvelU(j) * gamijk(k,i,j)
          enddo
       enddo
    enddo
    term2a = term2a + dtuUi
    !
    ! build the ugly long derivative term
    ! and get RUUsum:
    !  R_munu u^mu u^nu = R00 u^0^2 + 2 R0i u^0 u^i + Rij u^i u^j
    !
    !   --> and lower index of a^k, calculate magnitude
    !
    RUUsum    = fourR00 * (lorentz / alpijk)**2
    nabu_nabu = term1**2
    !
    ! zero-sum part of magnitude of accel
    maga = - (alpijk * aU0)**2
    do i=1,3
       RUUsum    = RUUsum + 2._dp * fourR0i(i) * fvelU(i) * lorentz / alpijk
       nabu_nabu = nabu_nabu + 2._dp * term2a(i) * term2b(i)

       do j=1,3
          RUUsum    = RUUsum + fourRij(i,j) * fvelU(i) * fvelU(j)
          nabu_nabu = nabu_nabu + term3a(i,j) * term3b(i,j)

          adi(i)    = adi(i) + gdown(i,j) * aUi(j)
       enddo
       maga = maga + aUi(i) * adi(i)
    enddo
    !maga = sqrt(maga)

    ! now get fluid R...
    fluidR = nabu_nabu - theta**2 + fourR + 2._dp * RUUsum

    ! store 4D components of aU and ad for output
    aUmu = (/ aU0, aUi(1), aUi(2), aUi(3) /)
    admu = (/ -alpijk**2*aU0, adi(1), adi(2), adi(3) /)

  end subroutine get_fluid_curvature




  !
  ! A subroutine to return the 4D Ricci tensor components
  !      --> from purely 3D, spatial objects
  !      --> to be called INSIDE a spatial loop
  !      --> ONLY called (from get_fluid_curvature) if tderivs .and. notearly
  !
  ! ** Also used for raytracer (called from interpolate.f90)
  !    --> in this case, we store our data differently. Current time is at START of t-dimension,
  !           whereas for regular mesc the current time is at END of t-dimension of arrs.
  !    --> SO change this so we pass in 'ti' = index of current time.
  !        - we assume this is either ti=1 or ti=nt, and will take time derivs of k accordingly
  !
  subroutine get_4D_ricci_tensor(nx,nt,ti,dx,dt,ipos,jpos,kpos,gdown_atstencil,gup,kij,&
       & kdown_atstencil,tracek,alp,threeRij,gamijk,dialp,kud_atstencil,R00,R0i,Rij,R)
    integer, intent(in) :: nx,nt,ti,ipos,jpos,kpos
    real(c_double), intent(in) :: dx,dt,dialp(3)
    real(c_double), intent(in) :: gup(3,3),threeRij(3,3),gamijk(3,3,3) ! 3R_ij and Christoffels
    real(c_double), intent(in) :: kij(6,nx,nx,nx,nt),tracek(nx,nx,nx,nt)
    real(c_double), dimension(13,3,3),   intent(in) :: gdown_atstencil,kdown_atstencil
    real(c_double), dimension(nx,nx,nx), intent(in) :: alp

    ! Output: the components of R_munu (symmetric) and 4D Ricci scalar
    real(c_double), intent(out) :: R00,R0i(3),Rij(3,3),R
    real(c_double), dimension(13,3,3), intent(out) :: kud_atstencil ! for violation_inloop

    real(c_double), dimension(3,3,nt) :: kdown_tstencil
    real(c_double), dimension(3,3) :: gupip1,gupim1,gupjp1,gupjm1,gupkp1,gupkm1
    real(c_double), dimension(3,3) :: gupip2,gupim2,gupjp2,gupjm2,gupkp2,gupkm2
    real(c_double), dimension(3,3) :: kdown,kud,kudip1,kudim1,kudjp1,kudjm1,kudkp1,kudkm1
    real(c_double), dimension(3,3) :: kudip2,kudim2,kudjp2,kudjm2,kudkp2,kudkm2

    real(c_double) :: didjalp(3,3),diK(3),digupij(3),djKudji(3)
    real(c_double) :: dxgup,dygup,dzgup,dxKud,dyKud,dzKud
    real(c_double) :: dtK,dtKij

    real(c_double) :: kUdKddsum(3,3)

    real(c_double) :: alpijk,alpip1,alpim1,alpjp1,alpjm1,alpkp1,alpkm1
    real(c_double) :: alpip2,alpim2,alpjp2,alpjm2,alpkp2,alpkm2
    real(c_double) :: trkijk,trkip1,trkim1,trkjp1,trkjm1,trkkp1,trkkm1
    real(c_double) :: trkip2,trkim2,trkjp2,trkjm2,trkkp2,trkkm2
    real(c_double) :: detgijk,ht
    !real(c_double) :: at,rhot,asq,rhostar,flrwhub,adot,addot,alpKUUKdd

    integer :: l,m,n
    integer :: ip1,jp1,kp1,im1,jm1,km1,ip2,im2,jp2,jm2,kp2,km2 ! for periodic boundaries
    integer :: ti_nm1,ti_nm2,ti_nm3,ti_nm4

    character(len=100) :: loc
    loc = "  get_4D_ricci_tensor"

    call apply_periodic(ipos,ip1,im1,nx)
    call apply_periodic(jpos,jp1,jm1,nx)
    call apply_periodic(kpos,kp1,km1,nx)
    call apply_periodic_fourth(ipos,ip2,im2,nx)
    call apply_periodic_fourth(jpos,jp2,jm2,nx)
    call apply_periodic_fourth(kpos,kp2,km2,nx)

    !
    ! We will pass in 'ti' being the index of the 'current' time in t-dim of arrays
    !    This is assumed to be either:
    !      1. ti = 1  --> raytracer, and take FORWARD derivatives (dt<0 in this case)
    !      2. ti = nt --> regular mesc, take BACKWARD derivatives (dt>0 in this case)
    ! below we set up the stencil for time derivatives
    ! in case 1. here we send to a BACKWARD deriv routine with dt --> -dt, (to avoid using case() statements too much in this routine)
    !select case(ti)
    !case(1)
    ! remove some compiler errors by initialising these
    ti_nm1 = -1; ti_nm2 = -1; ti_nm3 = -1; ti_nm4 = -1; ht = -1
    if (ti==1) then
        ! raytracer -- set up (f,fp1,fp2,fp3,fp4) and set ht = -dt so we can send to backward deriv routine
        ti_nm1 = 2
        ti_nm2 = 3
        ti_nm3 = 4
        ti_nm4 = 5
        ht     = -dt ! step for deriv routine
    !case(nt) <-- this doesn't work for some reason
    elseif (ti==nt) then
        ! regular mesc -- set up (f,fm1,fm2,fm3,fm4) and set ht = dt to send to backward deriv
        ti_nm1 = nt-1
        ti_nm2 = nt-2
        ti_nm3 = nt-3
        ti_nm4 = nt-4
        ht     = dt  ! setp for deriv routine
    !case default
    else
        call print_error("Please adjust your assumption that ti=1 or ti=nt, or fix your bug.",2,loc)
    endif ! select

    ! ------------------------------------------------------------------
    !
    ! Get spatial/time stencils we need for derivatives
    !
    ! ------------------------------------------------------------------
    !
    !
    ! Invert g_ij at all positions in stencil EXCEPT current (we have gup)
    !
    ! *_atstencil indices are:
    !     --> (ipos,jpos,kpos),im2,im1,ip1,ip2,jm2,jm1,jp1,jp2,km2,km1,kp1,kp2
    !
    call inv3x3(gdown_atstencil(2,:,:),gupim2,detgijk)
    call inv3x3(gdown_atstencil(5,:,:),gupip2,detgijk)
    call inv3x3(gdown_atstencil(6,:,:),gupjm2,detgijk)
    call inv3x3(gdown_atstencil(9,:,:),gupjp2,detgijk)
    call inv3x3(gdown_atstencil(10,:,:),gupkm2,detgijk)
    call inv3x3(gdown_atstencil(13,:,:),gupkp2,detgijk)
    call inv3x3(gdown_atstencil(3,:,:),gupim1,detgijk)
    call inv3x3(gdown_atstencil(4,:,:),gupip1,detgijk)
    call inv3x3(gdown_atstencil(7,:,:),gupjm1,detgijk)
    call inv3x3(gdown_atstencil(8,:,:),gupjp1,detgijk)
    call inv3x3(gdown_atstencil(11,:,:),gupkm1,detgijk)
    call inv3x3(gdown_atstencil(12,:,:),gupkp1,detgijk)

    !
    ! Calc K^i_j at all points in stencil
    !
    kdown = kdown_atstencil(1,:,:) ! K_ij at ipos,jpos,kpos
    call calc_up_down(kdown,gup,kud)
    call calc_up_down(kdown_atstencil(4,:,:),gupip1,kudip1)
    call calc_up_down(kdown_atstencil(3,:,:),gupim1,kudim1)
    call calc_up_down(kdown_atstencil(8,:,:),gupjp1,kudjp1)
    call calc_up_down(kdown_atstencil(7,:,:),gupjm1,kudjm1)
    call calc_up_down(kdown_atstencil(12,:,:),gupkp1,kudkp1)
    call calc_up_down(kdown_atstencil(11,:,:),gupkm1,kudkm1)
    ! Want fourth stencil anyway, since this is stored in _atstencil
    call calc_up_down(kdown_atstencil(5,:,:),gupip2,kudip2)
    call calc_up_down(kdown_atstencil(2,:,:),gupim2,kudim2)
    call calc_up_down(kdown_atstencil(9,:,:),gupjp2,kudjp2)
    call calc_up_down(kdown_atstencil(6,:,:),gupjm2,kudjm2)
    call calc_up_down(kdown_atstencil(13,:,:),gupkp2,kudkp2)
    call calc_up_down(kdown_atstencil(10,:,:),gupkm2,kudkm2)
    ! and output to be passed into violation calc!
    !      (can't think of a "nicer" way to code this other than every line explicitly...)
    kud_atstencil(1,:,:) = kud
    kud_atstencil(2,:,:) = kudim2
    kud_atstencil(3,:,:) = kudim1
    kud_atstencil(4,:,:) = kudip1
    kud_atstencil(5,:,:) = kudip2
    kud_atstencil(6,:,:) = kudjm2
    kud_atstencil(7,:,:) = kudjm1
    kud_atstencil(8,:,:) = kudjp1
    kud_atstencil(9,:,:) = kudjp2
    kud_atstencil(10,:,:) = kudkm2
    kud_atstencil(11,:,:) = kudkm1
    kud_atstencil(12,:,:) = kudkp1
    kud_atstencil(13,:,:) = kudkp2

    alpijk = alp(ipos,jpos,kpos)
    alpip1 = alp(ip1,jpos,kpos)
    alpim1 = alp(im1,jpos,kpos)
    alpjp1 = alp(ipos,jp1,kpos)
    alpjm1 = alp(ipos,jm1,kpos)
    alpkp1 = alp(ipos,jpos,kp1)
    alpkm1 = alp(ipos,jpos,km1)
    ! ti is the index of the 'current' time in t-dimension
    trkijk = tracek(ipos,jpos,kpos,ti)
    trkip1 = tracek(ip1,jpos,kpos,ti)
    trkim1 = tracek(im1,jpos,kpos,ti)
    trkjp1 = tracek(ipos,jp1,kpos,ti)
    trkjm1 = tracek(ipos,jm1,kpos,ti)
    trkkp1 = tracek(ipos,jpos,kp1,ti)
    trkkm1 = tracek(ipos,jpos,km1,ti)

    alpip2 = alp(ip2,jpos,kpos)
    alpim2 = alp(im2,jpos,kpos)
    alpjp2 = alp(ipos,jp2,kpos)
    alpjm2 = alp(ipos,jm2,kpos)
    alpkp2 = alp(ipos,jpos,kp2)
    alpkm2 = alp(ipos,jpos,km2)

    trkip2 = tracek(ip2,jpos,kpos,ti)
    trkim2 = tracek(im2,jpos,kpos,ti)
    trkjp2 = tracek(ipos,jp2,kpos,ti)
    trkjm2 = tracek(ipos,jm2,kpos,ti)
    trkkp2 = tracek(ipos,jpos,kp2,ti)
    trkkm2 = tracek(ipos,jpos,km2,ti)

    !
    ! Get time stencil for K_ij
    !
    do n=1,nt
       call get_metric_at_pos(ipos,jpos,kpos,nx,kij(:,:,:,:,n),kdown_tstencil(:,:,n))
    enddo

    ! ---------------------------------------------------------------------------
    !
    ! Take some derivatives we need
    !     -- since we're calculating them explicitly rather than in a loop
    !     -- this way is less lines than having a case for EVERY value of i and j in a loop
    !
    ! ---------------------------------------------------------------------------
    !
    ! diK -- which derivative we take depends on i (dialp passed in)
    !
    diK(1)   = deriv1fourth(trkip1,trkip2,trkim1,trkim2,dx) ! \pd_x K
    diK(2)   = deriv1fourth(trkjp1,trkjp2,trkjm1,trkjm2,dx) ! \pd_y K
    diK(3)   = deriv1fourth(trkkp1,trkkp2,trkkm1,trkkm2,dx) ! \pd_z K

    !
    ! get didj(alp) -- need every possible combination
    !
    didjalp(1,1) = deriv2fourth(alpip1,alpip2,alpijk,alpim1,alpim2,dx) ! d2x alp
    didjalp(2,2) = deriv2fourth(alpjp1,alpjp2,alpijk,alpjm1,alpjm2,dx) ! d2y alp
    didjalp(3,3) = deriv2fourth(alpkp1,alpkp2,alpijk,alpkm1,alpkm2,dx) ! d2z alp

    didjalp(1,2) = deriv2_mixfourth(ipos,jpos,nx,alp(:,:,kpos),dx) ! dxdy alp
    didjalp(1,3) = deriv2_mixfourth(ipos,kpos,nx,alp(:,jpos,:),dx) ! dxdz alp
    didjalp(2,3) = deriv2_mixfourth(jpos,kpos,nx,alp(ipos,:,:),dx) ! dydz alp

    ! apply symmetries
    didjalp(2,1) = didjalp(1,2) ! dydx alp = dxdy alp
    didjalp(3,1) = didjalp(1,3) ! dzdx alp = dxdz alp
    didjalp(3,2) = didjalp(2,3) ! dzdy alp = dydz alp

    !
    ! get d_0 K (depending on order set in options.f90)
    !
    select case(nord_dt)
    case(2)
       ! Use second order backward difference
       !dtK = deriv1_bckwrd2nd(trkijk,tracek(ipos,jpos,kpos,nt-1),tracek(ipos,jpos,kpos,nt-2),dt)
       dtK = deriv1_bckwrd2nd(trkijk,tracek(ipos,jpos,kpos,ti_nm1),tracek(ipos,jpos,kpos,ti_nm2),ht)
    case(3)
       ! Use third order backward difference
       !dtK = deriv1_bckwrd3rd(trkijk,tracek(ipos,jpos,kpos,nt-1),tracek(ipos,jpos,kpos,nt-2),&
        !    & tracek(ipos,jpos,kpos,nt-3),dt)
        dtK = deriv1_bckwrd3rd(trkijk,tracek(ipos,jpos,kpos,ti_nm1),tracek(ipos,jpos,kpos,ti_nm2),&
             & tracek(ipos,jpos,kpos,ti_nm3),ht)
    case(4)
       ! Use fourth order backward difference
       !dtK = deriv1_bckwrd4th(trkijk,tracek(ipos,jpos,kpos,nt-1),tracek(ipos,jpos,kpos,nt-2),&
        !    & tracek(ipos,jpos,kpos,nt-3),tracek(ipos,jpos,kpos,nt-4),dt)
        dtK = deriv1_bckwrd4th(trkijk,tracek(ipos,jpos,kpos,ti_nm1),tracek(ipos,jpos,kpos,ti_nm2),&
             & tracek(ipos,jpos,kpos,ti_nm3),tracek(ipos,jpos,kpos,ti_nm4),ht)
    case default
       call print_error("Only nord_dt = 2,3,4 implemented. Re-set and re-run.",2,loc)
    end select

    ! ------------------------------------------------------------------
    ! ------------------------------------------------------------------
    !
    ! Initialise calculation of R_ij and loop over components + sums
    !
    ! ------------------------------------------------------------------
    ! ------------------------------------------------------------------

    R00 = alpijk * dtK       ! zero-sum term
    R0i = 0._dp; Rij = 0._dp ! no terms outside of loops (zero-sum terms have components)
    R = 0._dp

    !alpKUUKdd = 0._dp

    ! -----------------------------------
    ! Loop breakdown:
    !    --> 'j' component loop for R_ij
    !    --> sum loop for R_00
    !    --> 'i' component loop for R_0i
    ! -----------------------------------
    do n=1,3

       ! -----------------------------------------------------
       !
       ! Spatial derivatives:
       !
       !  \pd_i g^ij = \pd_x g^xi  + \pd_y g^yi  + \pd_z g^zi
       ! \pd_j K^j_i = \pd_x K^x_i + \pd_y K^y_i + \pd_z K^z_i
       !
       dxgup = deriv1fourth(gupip1(1,n),gupip2(1,n),gupim1(1,n),gupim2(1,n),dx) ! \pd_x g^{xi}
       dygup = deriv1fourth(gupjp1(2,n),gupjp2(2,n),gupjm1(2,n),gupjm2(2,n),dx) ! \pd_y g^{yi}
       dzgup = deriv1fourth(gupkp1(3,n),gupkp2(3,n),gupkm1(3,n),gupkm2(3,n),dx) ! \pd_z g^{zi}

       dxKud = deriv1fourth(kudip1(1,n),kudip2(1,n),kudim1(1,n),kudim2(1,n),dx) ! \pd_x K^x_i
       dyKud = deriv1fourth(kudjp1(2,n),kudjp2(2,n),kudjm1(2,n),kudjm2(2,n),dx) ! \pd_y K^y_i
       dzKud = deriv1fourth(kudkp1(3,n),kudkp2(3,n),kudkm1(3,n),kudkm2(3,n),dx) ! \pd_z K^z_i
       ! sum each derivative term together
       digupij(n) = dxgup + dygup + dzgup
       djKudji(n) = dxKud + dyKud + dzKud

       ! -----------------------------------------------------

       R00 = R00 + alpijk * dialp(n) * digupij(n)

       ! ZERO-sum terms for R_0i
       R0i(n) = alpijk * (diK(n) - djKudji(n))

       ! ----------------------------
       ! 'i' component loop for R_ij
       ! sum loop for R_00, R_0i
       ! ----------------------------
       do m=1,3

          ! ------------------------------------------------
          !
          ! Time derivative: \pd_0 K_{m,n}
          !
          select case(nord_dt)
          case(2)
             ! Use second order backward difference
             !dtKij = deriv1_bckwrd2nd(kdown_tstencil(m,n,nt),kdown_tstencil(m,n,nt-1),kdown_tstencil(m,n,nt-2),dt)
             dtKij = deriv1_bckwrd2nd(kdown_tstencil(m,n,ti),kdown_tstencil(m,n,ti_nm1),kdown_tstencil(m,n,ti_nm2),ht)
          case(3)
             ! Use third order backward difference
             !dtKij = deriv1_bckwrd3rd(kdown_tstencil(m,n,nt),kdown_tstencil(m,n,nt-1),kdown_tstencil(m,n,nt-2),&
            !      & kdown_tstencil(m,n,nt-3),dt)
            dtKij = deriv1_bckwrd3rd(kdown_tstencil(m,n,ti),kdown_tstencil(m,n,ti_nm1),kdown_tstencil(m,n,ti_nm2),&
                 & kdown_tstencil(m,n,ti_nm3),ht)
          case(4)
             ! Use fourth order backward difference
             !dtKij = deriv1_bckwrd4th(kdown_tstencil(m,n,nt),kdown_tstencil(m,n,nt-1),kdown_tstencil(m,n,nt-2),&
            !      & kdown_tstencil(m,n,nt-3),kdown_tstencil(m,n,nt-4),dt)
            dtKij = deriv1_bckwrd4th(kdown_tstencil(m,n,ti),kdown_tstencil(m,n,ti_nm1),kdown_tstencil(m,n,ti_nm2),&
                 & kdown_tstencil(m,n,ti_nm3),kdown_tstencil(m,n,ti_nm4),ht)
          case default
             call print_error("Only nord_dt = 2,3,4 implemented. Re-set and re-run.",2,loc)
          end select
          ! ------------------------------------------------

          Rij(m,n) = threeRij(m,n) + trkijk * kdown(m,n) - (dtKij + didjalp(m,n)) / alpijk

          R00 = R00 + alpijk * gup(m,n) * didjalp(m,n) - alpijk**2 * kud(m,n) * kud(n,m)

          !alpKUUKdd = alpKUUKdd + alpijk**2 * kud(m,n) * kud(n,m)

          ! ------------------------------
          ! sum loop for R_ij, R_00, R_0i
          ! ------------------------------
          kUdKddsum(m,n) = 0._dp
          do l=1,3

             R00      = R00 + alpijk * gamijk(l,m,l) * gup(m,n) * dialp(n)

             R0i(n)   = R0i(n) + alpijk * ( gamijk(m,l,n) * kud(l,m) - gamijk(l,m,l) * kud(m,n) )

             Rij(m,n) = Rij(m,n) - 2._dp * kud(l,n) * kdown(m,l) + dialp(l) * gamijk(l,m,n) / alpijk

             kUdKddsum(m,n) = kUdKddsum(m,n) + 2._dp * kud(l,n) * kdown(m,l)

          enddo

          R = R + gup(m,n) * Rij(m,n)

       enddo
    enddo

    R = R - R00 / alpijk**2

  end subroutine get_4D_ricci_tensor



  !
  ! a function to return the fluid rest-frame spatial curvature
  ! as defined in Buchert+2018's new formalism -- using Hamiltonian-type like constraint
  ! \mathcal{R} = 16\pi G \rho - 2/3 \Theta^2 + 2 \sigma^2 - 2 w^2 (vorticity)
  !
  real(c_double) function fluid_restframe_curvature(rho,theta,sigma2,w2)
    real(c_double) :: theta,sigma2,w2,rho

    fluid_restframe_curvature = 16._dp * pi * rho - (2._dp / 3._dp) * theta**2 &
         + 2._dp * sigma2 - 2._dp * w2

  end function fluid_restframe_curvature





  !
  ! a subroutine to get the shear w.r.t the projection vector b_mu,nu = g_mu,nu + u_mu u_nu
  ! and to also get the fluid expansion \Theta (or \theta if tderivs=False)
  ! AND also get the vorticity w_munu (which should be zero for us; but we need to make sure!)
  !
  ! i.e. the shear and expansion in the new Buchert+ 2018/2019 formalism
  !
  ! * NOTE * this routine uses a linear approximation for \Theta \approx \theta and \sigma_new \approx \sigma_old (true
  !          to a factor of v^2 -- according to Pierre) if we have set tderivs=.False.
  !          OR if we are at iteration < nord_dt+1
  !
  ! double checked 11/3/19 OK. THOROUGHLY checked June 2019. It's fine.
  !
  subroutine get_expansion_shear_vort(ipos,jpos,kpos,nx,nt,dx,dt,notearly,gij_atstencil,&
       & kdown,Chrsijk,alp,dtalp,gup,vel1,vel2,vel3,tracek,fvelU_atstencil,umud,dialp,dtudt,&
       & dtudi,diudt,sigma2,theta,w2,lorentz,dwnsigmunu,upsigmunu,wddmunu)
    integer, intent(in) :: ipos,jpos,kpos,nx,nt
    real(c_double), intent(in) :: dx,dt,tracek,dtalp,alp(nx,nx,nx),Chrsijk(3,3,3) ! Christoffel symbols
    real(c_double), dimension(nx,nx,nx,nt), intent(in) :: vel1,vel2,vel3

    real(c_double), dimension(3,3),    intent(in) :: kdown,gup        ! since we already invert 3x3 in main ricci loop
    real(c_double), dimension(13,3,3), intent(in) :: gij_atstencil    ! the metric g_{ij} at entire (4th order) stencil
    logical, intent(in) :: notearly

    real(c_double), intent(out) :: dwnsigmunu(4,4),upsigmunu(4,4),sigma2,theta,w2,lorentz         ! w2 is vorticity scalar
    real(c_double), intent(out) :: fvelU_atstencil(13,3),umud(4),dialp(3),dtudt,dtudi(3),diudt(3) ! things passed into fluidR calc
    ! Need an additional variable here because if its not passed in, we don't want to set it... (using dwnwmunu(4,4) we did this)
    real(c_double), intent(out), optional :: wddmunu(4,4)!w0y,wxy
    !
    ! note: indices in names here follow the user guide (and written notes)
    !
    real(c_double) :: ud(3)!,upsigmunu(4,4)!,dwnsigmunu(4,4)
    real(c_double) :: upwmunu(4,4),dwnwmunu(4,4)
    real(c_double) :: gamsum1,gamsum2(3),gamsum3(3,3)   ! 1: \Gamma^d_{00} u_d, 2: \Gamma^d_{0l} u_d, 3: \Gamma^d_{lm} u_d
    real(c_double) :: bUdij(3,3),bUd0j(3),bddij         ! projection tensor b^i_j = \delta^i_j + u^i u_j, and b^0_j = -u^0 u_j, and b_ij
    real(c_double) :: bUd00,bUdi0(3),bdd00,bdd0i(3)
    real(c_double) :: deltaij,diudj(3,3)!,diudt(3)       ! \delta^i_j, \partial_i u_j, and \partial_i u_t
    real(c_double) :: gamjji_uUi_sum,dialp_uUi_sum      ! summation terms in \Theta
    real(c_double) :: gup_diudj_sum,gup_gamkij_udk_sum  ! summation terms in \theta
    real(c_double) :: gamjji,gamkij,gamilm,diuUi!,dtudt,dtudi(3)

    !
    ! things we need at stencil for spatial derivatives
    !
    real(c_double) :: alpijk,alpip1,alpim1,alpjp1,alpjm1,alpkp1,alpkm1
    real(c_double) :: alpip2,alpim2,alpjp2,alpjm2,alpkp2,alpkm2
    real(c_double) :: dxalp,dyalp,dzalp,dxuUx,dyuUy,dzuUz!,dialp(3)
    real(c_double) :: vsq,vsqip1,vsqim1,vsqjp1,vsqjm1,vsqkp1,vsqkm1
    real(c_double) :: vsqip2,vsqim2,vsqjp2,vsqjm2,vsqkp2,vsqkm2
    real(c_double) :: lorentzip1,lorentzim1,lorentzjp1,lorentzjm1,lorentzkp1,lorentzkm1
    real(c_double) :: lorentzip2,lorentzim2,lorentzjp2,lorentzjm2,lorentzkp2,lorentzkm2
    real(c_double), dimension(3) :: vels,velsip1,velsim1,velsjp1,velsjm1,velskp1,velskm1
    real(c_double), dimension(3) :: dtvelsi,velsnm1,velsnm2,velsnm3,velsnm4
    real(c_double), dimension(3) :: velsip2,velsim2,velsjp2,velsjm2,velskp2,velskm2
    real(c_double), dimension(3) :: uU,uUip1,uUim1,uUjp1,uUjm1,uUkp1,uUkm1
    real(c_double), dimension(3) :: uUip2,uUim2,uUjp2,uUjm2,uUkp2,uUkm2
    real(c_double), dimension(3) :: udip1,udim1,udjp1,udjm1,udkp1,udkm1!,ud
    real(c_double), dimension(3) :: udip2,udim2,udjp2,udjm2,udkp2,udkm2
    real(c_double) :: uUt,uUtip1,uUtim1,uUtjp1,uUtjm1,uUtkp1,uUtkm1
    real(c_double) :: uUtip2,uUtim2,uUtjp2,uUtjm2,uUtkp2,uUtkm2
    real(c_double) :: udt,udtip1,udtim1,udtjp1,udtjm1,udtkp1,udtkm1
    real(c_double) :: udtip2,udtim2,udtjp2,udtjm2,udtkp2,udtkm2,dtvsq

    real(c_double) :: wterm1,wterm2

    real(c_double), dimension(3,3) :: gdown,gdownip1,gdownim1,gdownjp1,gdownjm1,gdownkp1,gdownkm1
    real(c_double), dimension(3,3) :: gdownip2,gdownim2,gdownjp2,gdownjm2,gdownkp2,gdownkm2
    integer :: ip1,jp1,kp1,im1,jm1,km1,ip2,im2,jp2,jm2,kp2,km2 ! for periodic boundaries
    integer :: ii,jj,k,l,m,n                                   ! component and summation indices

    character(len=100) :: loc
    loc = "  get_expansion_shear_vort"

    ! --------------------------------------------------------
    !
    ! Get stencils for spatial deriv of alp and vels (for u_i)
    !
    ! --------------------------------------------------------

    call apply_periodic(ipos,ip1,im1,nx)
    call apply_periodic(jpos,jp1,jm1,nx)
    call apply_periodic(kpos,kp1,km1,nx)
    call apply_periodic_fourth(ipos,ip2,im2,nx)
    call apply_periodic_fourth(jpos,jp2,jm2,nx)
    call apply_periodic_fourth(kpos,kp2,km2,nx)
    vels = (/ vel1(ipos,jpos,kpos,nt), vel2(ipos,jpos,kpos,nt), vel3(ipos,jpos,kpos,nt) /)

    ! Time stencil - minimum 2nd order time derivs --> minimum 3 steps TOTAL (including t_n)
    if (tderivs .and. notearly) then
       if (nord_dt<2) call print_error(" Minimum required nord_dt = 2. Please re-set and re-run. ",2,loc)
       velsnm1 = (/ vel1(ipos,jpos,kpos,nt-1), vel2(ipos,jpos,kpos,nt-1), vel3(ipos,jpos,kpos,nt-1) /)
       velsnm2 = (/ vel1(ipos,jpos,kpos,nt-2), vel2(ipos,jpos,kpos,nt-2), vel3(ipos,jpos,kpos,nt-2) /)
       if (nord_dt >= 3) then
          !
          ! We want at least a 3rd order scheme, so need minimum 4 steps total (including t_n)
          !
          velsnm3 = (/ vel1(ipos,jpos,kpos,nt-3), vel2(ipos,jpos,kpos,nt-3), vel3(ipos,jpos,kpos,nt-3) /)
          if (nord_dt >= 4) then
             !
             ! We want at least a 4th order scheme, so need minumum 5 steps total (incluidng t_n)
             !
             velsnm4 = (/ vel1(ipos,jpos,kpos,nt-4), vel2(ipos,jpos,kpos,nt-4), vel3(ipos,jpos,kpos,nt-4) /)
             if (nord_dt > 4) call print_error(" nord_dt > 4 NOT IMPLEMENTED. Please set nord_dt <= 4 and re-run. ",2,loc)
          endif
       endif
    endif

    ! Spatial stencil
    alpijk = alp(ipos,jpos,kpos)
    alpip1 = alp(ip1,jpos,kpos)
    alpim1 = alp(im1,jpos,kpos)
    alpjp1 = alp(ipos,jp1,kpos)
    alpjm1 = alp(ipos,jm1,kpos)
    alpkp1 = alp(ipos,jpos,kp1)
    alpkm1 = alp(ipos,jpos,km1)
    velsip1 = (/ vel1(ip1,jpos,kpos,nt), vel2(ip1,jpos,kpos,nt), vel3(ip1,jpos,kpos,nt) /)
    velsim1 = (/ vel1(im1,jpos,kpos,nt), vel2(im1,jpos,kpos,nt), vel3(im1,jpos,kpos,nt) /)
    velsjp1 = (/ vel1(ipos,jp1,kpos,nt), vel2(ipos,jp1,kpos,nt), vel3(ipos,jp1,kpos,nt) /)
    velsjm1 = (/ vel1(ipos,jm1,kpos,nt), vel2(ipos,jm1,kpos,nt), vel3(ipos,jm1,kpos,nt) /)
    velskp1 = (/ vel1(ipos,jpos,kp1,nt), vel2(ipos,jpos,kp1,nt), vel3(ipos,jpos,kp1,nt) /)
    velskm1 = (/ vel1(ipos,jpos,km1,nt), vel2(ipos,jpos,km1,nt), vel3(ipos,jpos,km1,nt) /)

    alpip2 = alp(ip2,jpos,kpos)
    alpim2 = alp(im2,jpos,kpos)
    alpjp2 = alp(ipos,jp2,kpos)
    alpjm2 = alp(ipos,jm2,kpos)
    alpkp2 = alp(ipos,jpos,kp2)
    alpkm2 = alp(ipos,jpos,km2)
    velsip2 = (/ vel1(ip2,jpos,kpos,nt), vel2(ip2,jpos,kpos,nt), vel3(ip2,jpos,kpos,nt) /)
    velsim2 = (/ vel1(im2,jpos,kpos,nt), vel2(im2,jpos,kpos,nt), vel3(im2,jpos,kpos,nt) /)
    velsjp2 = (/ vel1(ipos,jp2,kpos,nt), vel2(ipos,jp2,kpos,nt), vel3(ipos,jp2,kpos,nt) /)
    velsjm2 = (/ vel1(ipos,jm2,kpos,nt), vel2(ipos,jm2,kpos,nt), vel3(ipos,jm2,kpos,nt) /)
    velskp2 = (/ vel1(ipos,jpos,kp2,nt), vel2(ipos,jpos,kp2,nt), vel3(ipos,jpos,kp2,nt) /)
    velskm2 = (/ vel1(ipos,jpos,km2,nt), vel2(ipos,jpos,km2,nt), vel3(ipos,jpos,km2,nt) /)

    ! gij_atstencil first dim is (pos,im2,im1,ip1,ip2,jm2,jm1,jp1,jp2,km2,km1,kp1,kp2)
    gdown = gij_atstencil(1,:,:)
    gdownim2 = gij_atstencil(2,:,:)
    gdownip2 = gij_atstencil(5,:,:)
    gdownjm2 = gij_atstencil(6,:,:)
    gdownjp2 = gij_atstencil(9,:,:)
    gdownkm2 = gij_atstencil(10,:,:)
    gdownkp2 = gij_atstencil(13,:,:)

    gdownim1 = gij_atstencil(3,:,:)
    gdownip1 = gij_atstencil(4,:,:)
    gdownjm1 = gij_atstencil(7,:,:)
    gdownjp1 = gij_atstencil(8,:,:)
    gdownkm1 = gij_atstencil(11,:,:)
    gdownkp1 = gij_atstencil(12,:,:)

    ! --------------------------------------------------------
    !
    ! Get lorentz factor at stencil
    !    1. Get v^2 at gij_atstencil
    !    2. Calc Lorentz factor at stencil from this
    !
    !    ALSO calculate dt(v^i) in here to save number of loops!
    !
    ! --------------------------------------------------------

    vsq = 0._dp ! \gamma_{ij}v^{i}v^{j} = v^2
    vsqip1=0._dp; vsqim1=0._dp; vsqjp1=0._dp
    vsqjm1=0._dp; vsqkp1=0._dp; vsqkm1=0._dp
    vsqip2=0._dp; vsqim2=0._dp; vsqjp2=0._dp
    vsqjm2=0._dp; vsqkp2=0._dp; vsqkm2=0._dp
    do k=1,3
       do l=1,3
          vsq = vsq + gdown(k,l) * vels(k) * vels(l)
          vsqip1 = vsqip1 + gdownip1(k,l) * velsip1(k) * velsip1(l)
          vsqim1 = vsqim1 + gdownim1(k,l) * velsim1(k) * velsim1(l)
          vsqjp1 = vsqjp1 + gdownjp1(k,l) * velsjp1(k) * velsjp1(l)
          vsqjm1 = vsqjm1 + gdownjm1(k,l) * velsjm1(k) * velsjm1(l)
          vsqkp1 = vsqkp1 + gdownkp1(k,l) * velskp1(k) * velskp1(l)
          vsqkm1 = vsqkm1 + gdownkm1(k,l) * velskm1(k) * velskm1(l)

          vsqip2 = vsqip2 + gdownip2(k,l) * velsip2(k) * velsip2(l)
          vsqim2 = vsqim2 + gdownim2(k,l) * velsim2(k) * velsim2(l)
          vsqjp2 = vsqjp2 + gdownjp2(k,l) * velsjp2(k) * velsjp2(l)
          vsqjm2 = vsqjm2 + gdownjm2(k,l) * velsjm2(k) * velsjm2(l)
          vsqkp2 = vsqkp2 + gdownkp2(k,l) * velskp2(k) * velskp2(l)
          vsqkm2 = vsqkm2 + gdownkm2(k,l) * velskm2(k) * velskm2(l)
       enddo

       ! --------------------------------------------------------
       !
       ! Take time derivs we need -- only done if we have tderivs=True
       !                     and we are on file number n >= nord_dt
       !
       ! --------------------------------------------------------
       if (tderivs .and. notearly) then
         ! First get dt(v^i) using data in memory
         select case(nord_dt)
         case(2)
            ! Use second order backward difference
            dtvelsi(k) = deriv1_bckwrd2nd(vels(k),velsnm1(k),velsnm2(k),dt)
         case(3)
            ! Use third order backward difference
            dtvelsi(k) = deriv1_bckwrd3rd(vels(k),velsnm1(k),velsnm2(k),velsnm3(k),dt)
         case(4)
            ! Use fourth order backward difference
            dtvelsi(k) = deriv1_bckwrd4th(vels(k),velsnm1(k),velsnm2(k),velsnm3(k),velsnm4(k),dt)
         case default
            call print_error("Only nord_dt = 2,3,4 implemented. Re-set and re-run.",2,loc)
         end select
       endif

     enddo
     lorentz     = 1._dp / sqrt( 1._dp - vsq )
     lorentzip1  = 1._dp / sqrt( 1._dp - vsqip1 )
     lorentzim1  = 1._dp / sqrt( 1._dp - vsqim1 )
     lorentzjp1  = 1._dp / sqrt( 1._dp - vsqjp1 )
     lorentzjm1  = 1._dp / sqrt( 1._dp - vsqjm1 )
     lorentzkp1  = 1._dp / sqrt( 1._dp - vsqkp1 )
     lorentzkm1  = 1._dp / sqrt( 1._dp - vsqkm1 )

     lorentzip2 = 1._dp / sqrt( 1._dp - vsqip2 )
     lorentzim2 = 1._dp / sqrt( 1._dp - vsqim2 )
     lorentzjp2 = 1._dp / sqrt( 1._dp - vsqjp2 )
     lorentzjm2 = 1._dp / sqrt( 1._dp - vsqjm2 )
     lorentzkp2 = 1._dp / sqrt( 1._dp - vsqkp2 )
     lorentzkm2 = 1._dp / sqrt( 1._dp - vsqkm2 )

    ! --------------------------------------------------------
    !
    ! get four-velocities u^mu and u_mu (by lowering index) at stencil
    !         (get u^i, u^0 and u_i, u_0 seperately)
    !
    ! --------------------------------------------------------

    uU    = lorentz * vels;       uUt = lorentz / alpijk ! u^i = \Gamma * v^i (from HydroBase defn), u^t = \Gamma / \alpha
    uUip1 = lorentzip1 * velsip1; uUtip1 = lorentzip1 / alpip1
    uUim1 = lorentzim1 * velsim1; uUtim1 = lorentzim1 / alpim1
    uUjp1 = lorentzjp1 * velsjp1; uUtjp1 = lorentzjp1 / alpjp1
    uUjm1 = lorentzjm1 * velsjm1; uUtjm1 = lorentzjm1 / alpjm1
    uUkp1 = lorentzkp1 * velskp1; uUtkp1 = lorentzkp1 / alpkp1
    uUkm1 = lorentzkm1 * velskm1; uUtkm1 = lorentzkm1 / alpkm1

    uUip2 = lorentzip2 * velsip2; uUtip2 = lorentzip2 / alpip2
    uUim2 = lorentzim2 * velsim2; uUtim2 = lorentzim2 / alpim2
    uUjp2 = lorentzjp2 * velsjp2; uUtjp2 = lorentzjp2 / alpjp2
    uUjm2 = lorentzjm2 * velsjm2; uUtjm2 = lorentzjm2 / alpjm2
    uUkp2 = lorentzkp2 * velskp2; uUtkp2 = lorentzkp2 / alpkp2
    uUkm2 = lorentzkm2 * velskm2; uUtkm2 = lorentzkm2 / alpkm2

    ! store u^i at stencil to be passed into fluidR routine for derivatives
    fvelU_atstencil(1,:) = uU
    fvelU_atstencil(2,:) = uUim2
    fvelU_atstencil(3,:) = uUim1
    fvelU_atstencil(4,:) = uUip1
    fvelU_atstencil(5,:) = uUip2
    fvelU_atstencil(6,:) = uUjm2
    fvelU_atstencil(7,:) = uUjm1
    fvelU_atstencil(8,:) = uUjp1
    fvelU_atstencil(9,:) = uUjp2
    fvelU_atstencil(10,:) = uUkm2
    fvelU_atstencil(11,:) = uUkm1
    fvelU_atstencil(12,:) = uUkp1
    fvelU_atstencil(13,:) = uUkp2

    !
    ! lower the indicies and calculate u_0 seperate
    !
    !  (AND calculate dt(v^2) to reduce loop numbers)
    !
    ud = 0._dp
    udip1 = 0._dp; udim1 = 0._dp; udjp1 = 0._dp
    udjm1 = 0._dp; udkp1 = 0._dp; udkm1 = 0._dp
    udip2 = 0._dp; udim2 = 0._dp; udjp2 = 0._dp
    udjm2 = 0._dp; udkp2 = 0._dp; udkm2 = 0._dp
    dtvsq = 0._dp
    do ii=1,3   ! component of u_i
       do l=1,3 ! sum loop
          ud(ii)    = ud(ii) + gdown(ii,l) * uU(l)
          udip1(ii) = udip1(ii) + gdownip1(ii,l) * uUip1(l)
          udim1(ii) = udim1(ii) + gdownim1(ii,l) * uUim1(l)
          udjp1(ii) = udjp1(ii) + gdownjp1(ii,l) * uUjp1(l)
          udjm1(ii) = udjm1(ii) + gdownjm1(ii,l) * uUjm1(l)
          udkp1(ii) = udkp1(ii) + gdownkp1(ii,l) * uUkp1(l)
          udkm1(ii) = udkm1(ii) + gdownkm1(ii,l) * uUkm1(l)

          udip2(ii) = udip2(ii) + gdownip2(ii,l) * uUip2(l)
          udim2(ii) = udim2(ii) + gdownim2(ii,l) * uUim2(l)
          udjp2(ii) = udjp2(ii) + gdownjp2(ii,l) * uUjp2(l)
          udjm2(ii) = udjm2(ii) + gdownjm2(ii,l) * uUjm2(l)
          udkp2(ii) = udkp2(ii) + gdownkp2(ii,l) * uUkp2(l)
          udkm2(ii) = udkm2(ii) + gdownkm2(ii,l) * uUkm2(l)

          if (tderivs .and. notearly) then
             ! Get dt(v^2)
             dtvsq = dtvsq + 2._dp * gdown(ii,l) * vels(ii) * dtvelsi(l) - &
                  2._dp * alpijk * kdown(ii,l) * vels(ii) * vels(l)
          endif
       enddo
    enddo

    udt    = - lorentz * alpijk  ! u_0 = - \Gamma \alpha
    udtip1 = - lorentzip1 * alpip1;    udtim1 = - lorentzim1 * alpim1
    udtjp1 = - lorentzjp1 * alpjp1;    udtjm1 = - lorentzjm1 * alpjm1
    udtkp1 = - lorentzkp1 * alpkp1;    udtkm1 = - lorentzkm1 * alpkm1

    udtip2 = - lorentzip2 * alpip2;    udtim2 = - lorentzim2 * alpim2
    udtjp2 = - lorentzjp2 * alpjp2;    udtjm2 = - lorentzjm2 * alpjm2
    udtkp2 = - lorentzkp2 * alpkp2;    udtkm2 = - lorentzkm2 * alpkm2

    !
    ! Store u_\mu for output
    umud = (/ udt, ud(1), ud(2), ud(3) /)

    ! --------------------------------------------------------
    !
    ! take spatial derivs of alp (3) and u^i (scalar)
    !
    ! --------------------------------------------------------

    ! \partial_i u^i
    dxuUx = deriv1fourth(uUip1(1),uUip2(1),uUim1(1),uUim2(1),dx)
    dyuUy = deriv1fourth(uUjp1(2),uUjp2(2),uUjm1(2),uUjm2(2),dx)
    dzuUz = deriv1fourth(uUkp1(3),uUkp2(3),uUkm1(3),uUkm2(3),dx)

    ! \partial \alp
    dxalp = deriv1fourth(alpip1,alpip2,alpim1,alpim2,dx)
    dyalp = deriv1fourth(alpjp1,alpjp2,alpjm1,alpjm2,dx)
    dzalp = deriv1fourth(alpkp1,alpkp2,alpkm1,alpkm2,dx)

    diuUi = dxuUx + dyuUy + dzuUz      ! \partial_i (u^i) -- sum
    dialp = (/ dxalp, dyalp, dzalp /)  ! \partial_i (alp)

    ! --------------------------------------------------------
    !
    ! build (3,3) tensors we need: \partial_i u_j, b^i_j
    !   build (3) vectors we need: b^0_j, b^i_0, b_0i, \partial_i u_t
    !       and b_munu components: b^0_0, b_00
    !
    ! and Christoffel sums: \Gamma^d_{00}u_d (scalar), \Gamma^d_{0l}u_d (3), \Gamma^d_{lm}u_d (3,3)
    ! and time derivatives: \partial_t u_t, \partial_t u_i below
    !
    ! --------------------------------------------------------

    ! initialise gamma sums with non-sum terms in them..
    !
    ! 1: \Gamma^d_{00}u_d = -\Gamma \partial_t(alp) + alp g^{ij} \partial_j(alp) u_i
    gamsum1 = - lorentz * dtalp

    ! 2: \Gamma^d_{0l}u_d = -\Gamma \partial_l(alp) - g^{ik} alp K_{kl} u_i
    gamsum2 = - lorentz * dialp ! vector (3)

    ! 3: \Gamma^d_{lm}u_d = \Gamma K_{lm} + \Gamma^i_{lm} u_i
    gamsum3 = lorentz * kdown   ! tensor (3,3)
    !
    ! b^0_0 = 1 - \Gam^2
    bUd00 = 1._dp - lorentz**2
    !
    ! b_00 = \alp^2 ( \Gamma^2 - 1 )
    bdd00 = alpijk**2 * ( lorentz**2 - 1._dp )
    !
    do jj=1,3
       !
       ! ** corrected - --> + May 8th, 2020
       ! b^0_j = u^0 u_i
       bUd0j(jj) = uUt * ud(jj)
       !
       ! b^i_0 = u^i u_0 = - \alp \Gamma u^i
       !bUdi0(jj) = - alpijk * lorentz * uU(jj)
       bUdi0(jj) = uU(jj) * udt
       !
       ! b_0i = u_0 u_i = - \alp^2 u^0 u_i
       !bdd0i(jj) = - alpijk**2 * uUt * ud(jj)
       bdd0i(jj) = udt * ud(jj)

       !
       ! \partial_i u_t
       select case(jj)
       case(1)
          ! x-derivative: \partial_x u_t
          diudt(jj) = deriv1fourth(udtip1,udtip2,udtim1,udtim2,dx)
       case(2)
          ! y-derivative: \partial_y u_t
          diudt(jj) = deriv1fourth(udtjp1,udtjp2,udtjm1,udtjm2,dx)
       case(3)
          ! z-derivative: \partial_z u_t
          diudt(jj) = deriv1fourth(udtkp1,udtkp2,udtkm1,udtkm2,dx)
       end select

       do ii=1,3

          ! Kronecker delta function
          if (ii==jj) then
             deltaij = 1._dp
          else
             deltaij = 0._dp
          endif

          ! b^i_j = d^i_j + u^i u_j
          bUdij(ii,jj) = deltaij + uU(ii) * ud(jj)

          ! \partial_ii u_jj
          select case(ii)
          case(1)
             ! x-derivative: \partial_x u_jj
             diudj(ii,jj) = deriv1fourth(udip1(jj),udip2(jj),udim1(jj),udim2(jj),dx)
          case(2)
             ! y-derivative: \partial_y u_jj
             diudj(ii,jj) = deriv1fourth(udjp1(jj),udjp2(jj),udjm1(jj),udjm2(jj),dx)
          case(3)
             ! z-derivative: \partial_z u_jj
             diudj(ii,jj) = deriv1fourth(udkp1(jj),udkp2(jj),udkm1(jj),udkm2(jj),dx)
          end select

          ! ------------------------------
          ! do \Gamma sums for \sigma_ij
          ! ------------------------------

          ! gamsum1 is a scalar and only has a double sum --> do it here
          gamsum1 = gamsum1 + alpijk * gup(ii,jj) * dialp(jj) * ud(ii)

          do l=1,3
             !
             ! gamsum3 has (3,3) components and a single sum --> need another loop
             gamilm = Chrsijk(l,ii,jj)
             gamsum3(ii,jj) = gamsum3(ii,jj) + gamilm * ud(l)

             !
             ! gamsum2 has (3) components and a double sum --> do it in here too
             gamsum2(jj) = gamsum2(jj) - gup(l,ii) * alpijk * kdown(ii,jj) * ud(l)

          enddo
       enddo
    enddo
    ! --------------------------------------------------------
    !
    ! sum loops outside component loops so we can calculate \Theta (or \theta)
    ! and the \sigma_ij terms with only one lot of sum loops
    !
    ! (also build Christoffel sums inside sum loops)
    !
    !    (AND get dtudt and dtudi to save loops!)
    !
    ! --------------------------------------------------------
    dtudt = - lorentz * dtalp ! we always have dtalp even if we don't have dt of v etc
    if (tderivs .and. notearly) then
      dtudt = dtudt + 0.5_dp * udt * lorentz**2 * dtvsq
    endif

    gamjji_uUi_sum = 0._dp; dialp_uUi_sum = 0._dp     ! \Theta sum terms
    gup_diudj_sum = 0._dp; gup_gamkij_udk_sum = 0._dp ! \theta sum terms
    !
    ! sum loops (component for dtudi)
    do l=1,3

       ! for \Theta
       ! \partial_i (alp) u^i / alp     !! single sum !!
       dialp_uUi_sum = dialp_uUi_sum + dialp(l) * uU(l) / alpijk

       ! (sum loop for dtudi)
       dtudi(l) = 0._dp
       do m=1,3

          ! for \Theta
          ! \Gamma^j_{ji} u^i           !! double sum !!
          gamjji = Chrsijk(l,l,m)
          gamjji_uUi_sum = gamjji_uUi_sum + gamjji * uU(m)

          ! for \theta
          ! g^{ij} \partial_i u_j       !! double sum !!
          gup_diudj_sum = gup_diudj_sum + gup(m,l) * diudj(m,l)

          do n=1,3
             ! for \theta
             ! g^{ij} \Gamma^k_{ij} u_k !! TRIPLE sum !!
             gamkij = Chrsijk(n,m,l)
             gup_gamkij_udk_sum = gup_gamkij_udk_sum + gup(m,l) * gamkij * ud(n)

          enddo

          if (tderivs .and. notearly) then
            ! Now get dtudi
            dtudi(l) = dtudi(l) + lorentz * gdown(l,m) * dtvelsi(m) + &
                 0.5_dp * lorentz**3 * gdown(l,m) * vels(m) * dtvsq - &
                 2._dp * alpijk * kdown(l,m) * uU(m)
          endif

       enddo

    enddo


    ! ----------------------------------------------------
    !
    ! Build expansion scalar \Theta (or \theta)
    !
    ! ----------------------------------------------------

    if (tderivs .and. notearly) then
       !
       ! we want \Theta (if we are on first few time steps, approximate with \theta below)
       !
       ! \Theta = \partial_i u^i - alp K u^t + \Gamma^j_{ji} u^i - u^t \partial_t(alp) / alp - \partial_t u_t / alp^2
       !             + 1/alp \partial_i(alp) u^i
       !
       theta = diuUi - alpijk * tracek * uUt + gamjji_uUi_sum - uUt * dtalp / alpijk - &
       & dtudt / alpijk**2 + dialp_uUi_sum

    else
       !
       ! we want \theta:
       !
       ! \theta = g^{ij} \partial_i u_j - K\Gamma - g^{ij} \Gamma^k_{ij} u_k
       !
       theta = gup_diudj_sum - tracek * lorentz - gup_gamkij_udk_sum

    endif

    ! ----------------------------------------------------
    !
    ! build shear \sigma_ij
    !
    ! NOTE: this loop will be too slow
    !       since calc_raised_comp also does 2x 1,3 loops
    !       incorporate this somehow into above loop?
    !
    ! ----------------------------------------------------

    !
    ! ---------------------
    ! calculate sigma(0,0) and w(0,0)
    !
    !     18 May, 2020:
    !     --> checked all gamsums and bUd,bdd's for typos vs. current mesc user guide
    !         (which was checked by Pierre, small typos fixed)
    !     --> sig_00, sig_0i, sig_ij forms are naively fine without corrections
    !
    ! ---------------------
    dwnwmunu = 0._dp; upwmunu = 0._dp
    dwnsigmunu = 0._dp; upsigmunu = 0._dp
    wterm1 = 0._dp; wterm2 = 0._dp
    if (tderivs .and. notearly) then
       !
       ! calculate sigma_00 using time derivs
       !        (w_00 = 0 TESTING)

       ! terms with no sums first
       dwnsigmunu(1,1) = bUd00**2 * (dtudt - gamsum1 ) - theta * bdd00 / 3._dp
       do l=1,3
          ! first (3) sum terms
          dwnsigmunu(1,1) = dwnsigmunu(1,1) + bUd00 * bUdi0(l) * ( dtudi(l) + diudt(l) - 2._dp * gamsum2(l) )
          do m=1,3
             ! next (3,3) sum terms
             dwnsigmunu(1,1) = dwnsigmunu(1,1) + 0.5_dp * bUdi0(l) * bUdi0(m) * ( diudj(l,m) + diudj(m,l) - 2._dp * gamsum3(l,m) )

             !===========================================================================================
             dwnwmunu(1,1) = dwnwmunu(1,1) + 0.5_dp * bUdi0(l) * bUdi0(m) * ( diudj(l,m) - diudj(m,l) )
             !===========================================================================================
          enddo
       enddo
       !
    else
       dwnsigmunu(1,1) = 0._dp
    endif

    !
    ! calc spatial parts of sig_munu and w_munu
    !
    do ii=1,3 ! component loop, sig,w components are ii+1
       !
       ! ------------------------------------------------------------
       ! calculate sigma(0,i) = sigma(i,0) and w(0,i) = - w(i,0)
       ! ------------------------------------------------------------

       if (tderivs .and. notearly) then
          !
          ! calculate sigma_0i, w_0i using time derivs
          !

          !
          ! terms with no sums first
          dwnsigmunu(1,ii+1) = bUd00 * bUd0j(ii) * (dtudt - gamsum1 ) - theta * bdd0i(ii) / 3._dp
          !
          do l=1,3
             ! then (3) sum terms
             dwnsigmunu(1,ii+1) = dwnsigmunu(1,ii+1) + 0.5_dp * ( dtudi(l) + diudt(l) - 2._dp * gamsum2(l) ) * &
                  ( bUd00 * bUdij(l,ii) + bUdi0(l) * bUd0j(ii) )

             ! w_0i
             dwnwmunu(1,ii+1)  = dwnwmunu(1,ii+1) + 0.5_dp * ( bUd00 * bUdij(l,ii) - bUdi0(l) * bUd0j(ii) ) * &
                  (dtudi(l) - diudt(l))

             do m=1,3
                ! next (3,3) sum term
                dwnsigmunu(1,ii+1) = dwnsigmunu(1,ii+1) + 0.5_dp * bUdi0(l) * bUdij(m,ii) * &
                     ( diudj(l,m) + diudj(m,l) - 2._dp * gamsum3(l,m) )

                dwnwmunu(1,ii+1)   = dwnwmunu(1,ii+1) + 0.5_dp * bUdi0(l) * bUdij(m,ii) * ( diudj(l,m) - diudj(m,l) )
             enddo
          enddo
          !
       else
          dwnsigmunu(1,ii+1) = 0._dp
          dwnwmunu(1,ii+1)   = 0._dp
       endif
       !
       ! symmetric: sigma_0i = sigma_i0
       ! anti-symmetric w_0i = - w_i0
       !
       dwnsigmunu(ii+1,1) = dwnsigmunu(1,ii+1)
       dwnwmunu(ii+1,1)   = - dwnwmunu(1,ii+1)
       !
       do jj=1,3
          ! ------------------------------------------------------------
          ! calculate sigma(i,j) and w(i,j)
          ! ------------------------------------------------------------
          !
          ! projection tensor
          bddij = gdown(ii,jj) + ud(ii) * ud(jj)
          !
          if (tderivs .and. notearly) then
             !
             ! calculate sigma_ij and w_ij using time derivs
             !

             ! no additional sums in this term
             dwnsigmunu(ii+1,jj+1) = bUd0j(ii) * bUd0j(jj) * ( dtudt - gamsum1 ) - theta * bddij / 3._dp
             !
             do l=1,3
                !
                ! single sum per component in this term
                dwnsigmunu(ii+1,jj+1) = dwnsigmunu(ii+1,jj+1) + 0.5_dp * ( diudt(l) + dtudi(l) - 2._dp * gamsum2(l) ) * &
                     ( bUd0j(ii) * bUdij(l,jj) + bUdij(l,ii) * bUd0j(jj) )

                dwnwmunu(ii+1,jj+1)   = dwnwmunu(ii+1,jj+1) + 0.5_dp * (bUd0j(ii) * bUdij(l,jj) - bUdij(l,ii) * bUd0j(jj)) * &
                    (dtudi(l) - diudt(l))

                if (ii==1 .and. jj==2) wterm1 = wterm1 + 0.5_dp * (bUd0j(ii) * bUdij(l,jj) - bUdij(l,ii) * bUd0j(jj)) * &
                    (dtudi(l) - diudt(l))
                !if (ii==1 .and. jj==2) print*, (dtudi(l) - diudt(l))
                !if (ii==1 .and. jj==2) print*, (bUd0j(ii) * bUdij(l,jj) - bUdij(l,ii) * bUd0j(jj))

                do m=1,3
                   !
                   ! double sum per component in this term
                   dwnsigmunu(ii+1,jj+1) = dwnsigmunu(ii+1,jj+1) + 0.5_dp * bUdij(l,ii) * bUdij(m,jj) * &
                        (diudj(l,m) + diudj(m,l) - 2._dp * gamsum3(l,m) )

                   dwnwmunu(ii+1,jj+1)   = dwnwmunu(ii+1,jj+1) + 0.5_dp * bUdij(l,ii) * bUdij(m,jj) * (diudj(l,m) - diudj(m,l))

                   ! get x,y component of this term
                   if (ii==1 .and. jj==2) wterm2 = wterm2 + 0.5_dp * bUdij(l,ii) * bUdij(m,jj) * (diudj(l,m) - diudj(m,l))


                enddo
             enddo
          else
             !
             ! we aren't running on all snapshots --> approx with "old" \sigma (from h_ij)
             !
             dwnsigmunu(ii+1,jj+1) = 0.5_dp * (diudj(ii,jj) + diudj(jj,ii) - 2._dp * gamsum3(ii,jj) ) - theta * gdown(ii,jj) / 3._dp
             dwnwmunu(ii+1,jj+1)   = 0._dp ! this is zero to O(3) (According to Pierre)
          endif
       enddo
    enddo
    !
    ! if testing, output omega_mu,nu
    !
    if (present(wddmunu)) wddmunu = dwnwmunu
    !if (present(w0y))     w0y = dwnwmunu(1,3)
    !if (present(wxy))     wxy = dwnwmunu(2,3)

    !
    ! raise indices of sigma_ij to calculate sigma^2
    ! --> we need another loop here because we need the ENTIRE (3,3) dwnsigij to do this
    ! -- future: calculate upsigmunu alongside dwnsigmunu (but would still need another loop to calc sigma^2 ... )
    !
    sigma2 = 0._dp; w2 = 0._dp
    do ii=1,4
       do jj=1,4

          !
          ! call calc_raised_comp_4D(i,j,Add,upgij,upg00,upg0i,Auu)
          call calc_raised_comp_4D(ii,jj,dwnsigmunu,gup,-alpijk**(-2),(/0._dp,0._dp,0._dp/),upsigmunu(ii,jj))
          call calc_raised_comp_4D(ii,jj,dwnwmunu,gup,-alpijk**(-2),(/0._dp,0._dp,0._dp/),upwmunu(ii,jj))
          !
          ! \sigma^2 = 1/2 sigma^{ij} sigma_{ij}
          sigma2 = sigma2 + 0.5_dp * upsigmunu(ii,jj) * dwnsigmunu(ii,jj)
          w2     = w2 + 0.5_dp * upwmunu(ii,jj) * dwnwmunu(ii,jj)

       enddo
    enddo

  end subroutine get_expansion_shear_vort


  !
  ! A cute little subroutine to get dialp in a spatial loop
  !       -- i.e. when we want to hide the derivative calls + the periodic bc implementation
  !
  subroutine get_dialp(nx,dx,ipos,jpos,kpos,alp,dialp)
      integer, intent(in) :: nx,ipos,jpos,kpos ! grid size and index position
      real(c_double), intent(in) :: dx,alp(nx,nx,nx)
      real(c_double), intent(out) :: dialp(3)

      real(c_double) :: dxalp,dyalp,dzalp
      real(c_double) :: alpip1,alpip2,alpim1,alpim2,alpjp1,&
        & alpjp2,alpjm1,alpjm2,alpkp1,alpkp2,alpkm1,alpkm2
      integer :: ip1,ip2,im1,im2,jp1,jp2,jm1,jm2,kp1,kp2,km1,km2

      call apply_periodic(ipos,ip1,im1,nx)
      call apply_periodic(jpos,jp1,jm1,nx)
      call apply_periodic(kpos,kp1,km1,nx)
      call apply_periodic_fourth(ipos,ip2,im2,nx)
      call apply_periodic_fourth(jpos,jp2,jm2,nx)
      call apply_periodic_fourth(kpos,kp2,km2,nx)

      alpip1 = alp(ip1,jpos,kpos)
      alpim1 = alp(im1,jpos,kpos)
      alpjp1 = alp(ipos,jp1,kpos)
      alpjm1 = alp(ipos,jm1,kpos)
      alpkp1 = alp(ipos,jpos,kp1)
      alpkm1 = alp(ipos,jpos,km1)

      alpip2 = alp(ip2,jpos,kpos)
      alpim2 = alp(im2,jpos,kpos)
      alpjp2 = alp(ipos,jp2,kpos)
      alpjm2 = alp(ipos,jm2,kpos)
      alpkp2 = alp(ipos,jpos,kp2)
      alpkm2 = alp(ipos,jpos,km2)

      !
      ! get \partial \alp
      !
      dxalp = deriv1fourth(alpip1,alpip2,alpim1,alpim2,dx)
      dyalp = deriv1fourth(alpjp1,alpjp2,alpjm1,alpjm2,dx)
      dzalp = deriv1fourth(alpkp1,alpkp2,alpkm1,alpkm2,dx)
      dialp = (/ dxalp, dyalp, dzalp /)  ! \partial_i (alp)

  end subroutine get_dialp



end module roots
