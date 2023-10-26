module riemann
    !
    ! A module to calculate components of the Riemann tensor
    !
    use options, only:c_double,nord_dt,clen
    use periodic, only:apply_periodic,apply_periodic_fourth
    use derivatives, only:deriv2_mixfourth,deriv2fourth,deriv1fourth,deriv1,&
        & deriv1_forward4th,deriv1_forward3rd,deriv1_forward2nd,deriv2_mix,deriv2,&
        & deriv1_bckwrd4th,deriv1_bckwrd3rd,deriv1_bckwrd2nd
    !use interpolate, only:nord_interp_t
    use manipulations, only:get_metric_at_pos
    use prints, only:print_error
    implicit none

contains

    !
    ! A subroutine to return the Weyl optical scalar \mathscr{W} = - 0.5 R_{AB} X^A X^B
    !
    !    ** checked HMac April 28th, 2021, OK
    !
    subroutine get_WWeyl(RAB,kU0,kUi,smu1,smu2,ImWWeyl,ReWWeyl)
        real(c_double), intent(in) :: RAB(6,6),kU0,kUi(3),smu1(4),smu2(4)

        real(c_double), intent(out) :: ImWWeyl,ReWWeyl ! Real and Imaginary parts of \mathscr{W}

        real(c_double) :: X1(6),X2(6) ! X_1^A \equiv k^mu s^nu_1        (similarly for X2)
        real(c_double) :: Y1(6),Y2(6) ! Y_1^B \equiv k^\lambda s^\rho_1 (similarly for Y2) -- note different index order
        integer :: A,B ! \in [01, 02, 03, 12, 13, 23]

        ! How to loop?
        !   we have, e.g. RAB(1,:) = 0101,0102,0103,0112,0113,0123
        !  so [A,B] \in (01,02,03,12,13,23)
        !  --> So, as long as we build X1, X2 with same ordering as we did RAB,
        !       we should just be able to loop over A,B..

        ! 1.
        ! Build X1 = X_1^[A] \equiv s^mu_1 k^nu - s^nu_1 k^mu
        ! * note we've left out the 1/2 in front of antisymmetric X,Y since sum has 2 outside R_AB
        ! (mu,nu) \in 01,02,03,12,13,23

        !
        ! ... seems this is easiest to just write out
        X1(1) = smu1(1) * kUi(1) - smu1(2) * kU0     ! 01: s^0_1 k^1 - s^1_1 k^0 (smu(4) so i+1)
        X1(2) = smu1(1) * kUi(2) - smu1(3) * kU0     ! 02: s^0_1 k^2 - s^2_1 k^0
        X1(3) = smu1(1) * kUi(3) - smu1(4) * kU0     ! 03: s^0_1 k^3 - s^3_1 k^0
        X1(4) = smu1(2) * kUi(2) - smu1(3) * kUi(1)  ! 12: s^1_1 k^2 - s^2_1 k^1
        X1(5) = smu1(2) * kUi(3) - smu1(4) * kUi(1)  ! 13: s^1_1 k^3 - s^3_1 k^1
        X1(6) = smu1(3) * kUi(3) - smu1(4) * kUi(2)  ! 23: s^2_1 k^3 - s^3_1 k^2
        !
        ! 2.
        ! Build X2 (same as X1 with s^mu_1 --> s^mu_2)
        X2(1) = smu2(1) * kUi(1) - smu2(2) * kU0     ! 01: s^0_2 k^1 - s^1_2 k^0
        X2(2) = smu2(1) * kUi(2) - smu2(3) * kU0     ! 02: s^0_2 k^2 - s^2_2 k^0
        X2(3) = smu2(1) * kUi(3) - smu2(4) * kU0     ! 03: s^0_2 k^3 - s^3_2 k^0
        X2(4) = smu2(2) * kUi(2) - smu2(3) * kUi(1)  ! 12: s^1_2 k^2 - s^2_2 k^1
        X2(5) = smu2(2) * kUi(3) - smu2(4) * kUi(1)  ! 13: s^1_2 k^3 - s^3_2 k^1
        X2(6) = smu2(3) * kUi(3) - smu2(4) * kUi(2)  ! 23: s^2_2 k^3 - s^3_2 k^2
        !
        ! 3.
        ! Build Y1 = Y_1^B \equiv k^lambda s^rho_1
        ! (lambda,rho) \in 01,02,03,12,13,23
        Y1(1) = kU0 * smu1(2) - kUi(1) * smu1(1)    ! 01: k^0 s^1_1 - k^1 s^0_1
        Y1(2) = kU0 * smu1(3) - kUi(2) * smu1(1)    ! 02: k^0 s^2_1 - k^2 s^0_1
        Y1(3) = kU0 * smu1(4) - kUi(3) * smu1(1)    ! 03: k^0 s^3_1 - k^3 s^0_1
        Y1(4) = kUi(1) * smu1(3) - kUi(2) * smu1(2) ! 12: k^1 s^2_1 - k^2 s^1_1
        Y1(5) = kUi(1) * smu1(4) - kUi(3) * smu1(2) ! 13: k^1 s^3_1 - k^3 s^1_1
        Y1(6) = kUi(2) * smu1(4) - kUi(3) * smu1(3) ! 23: k^2 s^3_1 - k^3 s^2_1
        !
        ! 4.
        ! Build Y2 (same as Y1 with s^mu_1 --> s^mu_2)
        Y2(1) = kU0 * smu2(2) - kUi(1) * smu2(1)    ! 01: k^0 s^1_2 - k^1 s^0_2
        Y2(2) = kU0 * smu2(3) - kUi(2) * smu2(1)    ! 02: k^0 s^2_2 - k^2 s^0_2
        Y2(3) = kU0 * smu2(4) - kUi(3) * smu2(1)    ! 03: k^0 s^3_2 - k^3 s^0_2
        Y2(4) = kUi(1) * smu2(3) - kUi(2) * smu2(2) ! 12: k^1 s^2_2 - k^2 s^1_2
        Y2(5) = kUi(1) * smu2(4) - kUi(3) * smu2(2) ! 13: k^1 s^3_2 - k^3 s^1_2
        Y2(6) = kUi(2) * smu2(4) - kUi(3) * smu2(3) ! 23: k^2 s^3_2 - k^3 s^2_2

        ImWWeyl = 0.d0; ReWWeyl = 0.d0
        do A=1,6
            do B=1,6

                ReWWeyl = ReWWeyl - RAB(A,B) * (X1(A) * Y1(B) - X2(A) * Y2(B))

                ImWWeyl = ImWWeyl + RAB(A,B) * (X1(A) * Y2(B) + X2(A) * Y1(B))
            enddo
        enddo


    end subroutine get_WWeyl



    !
    ! A subroutine to get R_{AB} = R_{mu,nu,lambda,rho} (Riemann tensor) for combinations
    !      of [A,B] we need to calculate the Weyl scalar \mathscr{W}
    !
    !   --> this is intended to be called from inside the interpolation loop
    !
    !  ** checked April 12th, 2021, OK HMAC
    !
    subroutine get_RAB(nx,nt,ti,ipos,jpos,kpos,dx,dt,gij,kij,alp,gamijk,gdown,gup,kdown,RAB)
        integer, intent(in) :: nx,nt,ti,ipos,jpos,kpos
        real(c_double), intent(in) :: dx,dt,alp(nx,nx,nx),gij(6,nx,nx,nx)
        real(c_double), intent(in) :: kij(6,nx,nx,nx,nt)
        real(c_double), intent(in) :: gamijk(3,3,3,nx,nx,nx)
        real(c_double), intent(in) :: gdown(3,3),kdown(3,3),gup(3,3)

        real(c_double), intent(out) :: RAB(6,6)

        ! -------------------------------------------
        !           Indices that we need:
        !    A,B \in [01, 02, 03, 12, 13, 23]
        ! -------------------------------------------

        !
        ! R0i0j terms are: [AB] = [0101,0102,0103,0202,0203,0303]
        !
        ! Call to get_R0i0j gives us the whole (3,3) for all i,j
        !    this is the top left corner of R_AB: i.e. R_{AB}(1:3,1:3)
        !     (symmetries are applied in  this routine)
        !

        call get_R0i0j(nx,nt,ti,dx,dt,ipos,jpos,kpos,kij,alp,gamijk,kdown,gup,RAB(1:3,1:3))

        !
        ! R_{0ijk} terms are: [AB] = [0112,0113,0123,0212,0213,0223,0312,0313,0323]
        !
        ! Get each of these individually and store in the upper right part of R_{AB}
        !      (couldn't think of a nice way to loop over these...
        !                  easier to do it explicitly for now)
        !
        ! R_{0112,0113,0123} = RAB(1,4:6)
        call get_R0ijk(nx,ipos,jpos,kpos,dx,kij(:,:,:,:,ti),alp,gamijk,kdown,1,1,2,RAB(1,4))
        call get_R0ijk(nx,ipos,jpos,kpos,dx,kij(:,:,:,:,ti),alp,gamijk,kdown,1,1,3,RAB(1,5))
        call get_R0ijk(nx,ipos,jpos,kpos,dx,kij(:,:,:,:,ti),alp,gamijk,kdown,1,2,3,RAB(1,6))
        ! Apply symmetries
        RAB(4,1) = RAB(1,4)
        RAB(5,1) = RAB(1,5)
        RAB(6,1) = RAB(1,6)
        !
        ! R_{0212,0213,0223} = RAB(2,4:6)
        call get_R0ijk(nx,ipos,jpos,kpos,dx,kij(:,:,:,:,ti),alp,gamijk,kdown,2,1,2,RAB(2,4))
        call get_R0ijk(nx,ipos,jpos,kpos,dx,kij(:,:,:,:,ti),alp,gamijk,kdown,2,1,3,RAB(2,5))
        call get_R0ijk(nx,ipos,jpos,kpos,dx,kij(:,:,:,:,ti),alp,gamijk,kdown,2,2,3,RAB(2,6))
        ! Apply symmetries
        RAB(4,2) = RAB(2,4)
        RAB(5,2) = RAB(2,5)
        RAB(6,2) = RAB(2,6)
        !
        ! R_{0312,0313,0323} = RAB(3,4:6)
        call get_R0ijk(nx,ipos,jpos,kpos,dx,kij(:,:,:,:,ti),alp,gamijk,kdown,3,1,2,RAB(3,4))
        call get_R0ijk(nx,ipos,jpos,kpos,dx,kij(:,:,:,:,ti),alp,gamijk,kdown,3,1,3,RAB(3,5))
        call get_R0ijk(nx,ipos,jpos,kpos,dx,kij(:,:,:,:,ti),alp,gamijk,kdown,3,2,3,RAB(3,6))
        ! Apply symmetries
        RAB(4,3) = RAB(3,4)
        RAB(5,3) = RAB(3,5)
        RAB(6,3) = RAB(3,6)
        !
        ! R_{ijkl} terms are: [AB] = [1212,1213,1223,1313,1323,2323]
        !
        ! Again, just get these explicitly...
        !
        ! R_{1212,1213,1223} = RAB(4,4:6)
        call get_Rijkl(nx,dx,ipos,jpos,kpos,gamijk,gdown,kdown,1,2,1,2,RAB(4,4))
        call get_Rijkl(nx,dx,ipos,jpos,kpos,gamijk,gdown,kdown,1,2,1,3,RAB(4,5))
        call get_Rijkl(nx,dx,ipos,jpos,kpos,gamijk,gdown,kdown,1,2,2,3,RAB(4,6))
        ! Apply symmetries
        RAB(5,4) = RAB(4,5)
        RAB(6,4) = RAB(4,6)
        !
        ! R_{1313,1323} = RAB(5,5:6)
        call get_Rijkl(nx,dx,ipos,jpos,kpos,gamijk,gdown,kdown,1,3,1,3,RAB(5,5))
        call get_Rijkl(nx,dx,ipos,jpos,kpos,gamijk,gdown,kdown,1,3,2,3,RAB(5,6))
        ! Apply symmetries
        RAB(6,5) = RAB(5,6)
        !
        ! R_{2323} = RAB(6,6)
        call get_Rijkl(nx,dx,ipos,jpos,kpos,gamijk,gdown,kdown,2,3,2,3,RAB(6,6))


    end subroutine get_RAB



    !
    ! Return R_{0i0j} at a particular position for all i,j
    !      (note this is symmetric s.t. only 6 indep components)
    !
    !    We did it this way rather than passing in indices i,j because for raytracer
    !        we need 6/9 components, so it's worth outputting the whole (3,3)
    !
    ! ** checked Mar 31 2021 - HMac OK
    !
    subroutine get_R0i0j(nx,nt,ti,dx,dt,ipos,jpos,kpos,kij,alp,gamijk,kdown,gup,R0i0j)
        integer, intent(in) :: nx,nt,ti,ipos,jpos,kpos ! spatial position
        real(c_double), intent(in) :: dx,dt,kij(6,nx,nx,nx,nt)
        real(c_double), intent(in) :: gamijk(3,3,3,nx,nx,nx)
        real(c_double), intent(in) :: alp(nx,nx,nx),kdown(3,3),gup(3,3)

        real(c_double), intent(out) :: R0i0j(3,3)

        real(c_double), dimension(3,3) :: didjalp
        real(c_double) :: alpijk,dxalp,dyalp,dzalp,dialp(3)
        real(c_double) :: alpjp1,alpjp2,alpjm1,alpjm2,alpim1,alpim2
        real(c_double) :: alpkp1,alpkp2,alpkm1,alpkm2,alpip1,alpip2
        real(c_double) :: dkalp_gamsum,KKsum,Gamma,dtKij,dtf

        integer :: ip1,im1,jp1,jm1,kp1,km1,ip2,im2,jp2,jm2,kp2,km2
        integer :: tip1,tip2,tip3,tip4
        integer :: ii,jj,i,k,l,ivals(6),jvals(6)
        character(len=clen), parameter :: loc = " get_R0i0j"

        ! First; sort out time indices for taking forward time derivative
        !        of Kij. This routine is either called from RT or "regular" mesc
        !        these both have different ordering when storing data in the time dimensino of arrays
        !        (this should be fixed to be consistent in future... maybe)
        ! Where our call is coming from will be encoded in 'ti' (index of "current time" in time dim)
        ! 1. RT:   ti = 1, current time at beginning of arrays
        ! 2. mesc: ti = nt, current time at end of arrays
        !
        ! which case we have will depend how we take time deriv of Kij. We set this up as follows:
        if (ti==1) then
            ! We are being called from the raytracer, these are straightforward
            ! ti=1 is current time, ti=2 is t-dt, ti=3 is t-2dt (backwards ordered)
            !   BUT we have already set dt-->-dt, so to mesc this is "forward" in time
            tip1 = ti + 1
            tip2 = ti + 2
            tip3 = ti + 3
            tip4 = ti + 4
            dtf  = dt
        elseif (ti==nt) then
            ! We are being called from regular mescaline (via ricci.f90)
            !    these are backwards w.r.t RT
            ! ti=nt is current time, nt-1 is t-dt, nt-2 is t-2dt,... (forwards ordered)
            ! set dt-->-dt so we can use forward difference for dt Kij
            tip1 = ti - 1
            tip2 = ti - 2
            tip3 = ti - 3
            tip4 = ti - 4
            dtf  = -dt
        else
            ! save some compilation errors
            tip1 = -1; tip2 = -1; tip3 = -1; tip4 = -1; dtf = -1
            call print_error("Unrecognised ordering of time dimension. Check yo' self.",2,loc)
        endif

        !
        ! Apply periodic boundary conditions
        call apply_periodic(ipos,ip1,im1,nx)
        call apply_periodic(jpos,jp1,jm1,nx)
        call apply_periodic(kpos,kp1,km1,nx)
        call apply_periodic_fourth(ipos,ip2,im2,nx)
        call apply_periodic_fourth(jpos,jp2,jm2,nx)
        call apply_periodic_fourth(kpos,kp2,km2,nx)

        !
        ! Get lapse at spatial stencil and take derivatives
        !
        alpijk = alp(ipos,jpos,kpos)
        alpip1 = alp(ip1,jpos,kpos); alpim1 = alp(im1,jpos,kpos)
        alpjp1 = alp(ipos,jp1,kpos); alpjm1 = alp(ipos,jm1,kpos)
        alpkp1 = alp(ipos,jpos,kp1); alpkm1 = alp(ipos,jpos,km1)
        alpip2 = alp(ip2,jpos,kpos); alpim2 = alp(im2,jpos,kpos)
        alpjp2 = alp(ipos,jp2,kpos); alpjm2 = alp(ipos,jm2,kpos)
        alpkp2 = alp(ipos,jpos,kp2); alpkm2 = alp(ipos,jpos,km2)

        didjalp(1,1) = deriv2fourth(alpip1,alpip2,alpijk,alpim1,alpim2,dx) ! d2x alp
        didjalp(2,2) = deriv2fourth(alpjp1,alpjp2,alpijk,alpjm1,alpjm2,dx) ! d2y alp
        didjalp(3,3) = deriv2fourth(alpkp1,alpkp2,alpijk,alpkm1,alpkm2,dx) ! d2z alp

        didjalp(1,2) = deriv2_mixfourth(ipos,jpos,nx,alp(:,:,kpos),dx) ! dxdy alp
        didjalp(1,3) = deriv2_mixfourth(ipos,kpos,nx,alp(:,jpos,:),dx) ! dxdz alp
        didjalp(2,3) = deriv2_mixfourth(jpos,kpos,nx,alp(ipos,:,:),dx) ! dydz alp

        ! \partial \alp
        dxalp = deriv1fourth(alpip1,alpip2,alpim1,alpim2,dx)
        dyalp = deriv1fourth(alpjp1,alpjp2,alpjm1,alpjm2,dx)
        dzalp = deriv1fourth(alpkp1,alpkp2,alpkm1,alpkm2,dx)

        ! apply symmetries
        didjalp(2,1) = didjalp(1,2) ! dydx alp = dxdy alp
        didjalp(3,1) = didjalp(1,3) ! dzdx alp = dxdz alp
        didjalp(3,2) = didjalp(2,3) ! dzdy alp = dydz alp
        dialp = (/ dxalp, dyalp, dzalp /)

       !
       ! Loop like this only x6 rather than 3x3=9 every time its called, apply symmetries after
       !
       ivals = (/ 1, 1, 1, 2, 2, 3 /)
       jvals = (/ 1, 2, 3, 2, 3, 3 /)

       R0i0j = 0.d0
       do i=1,6
           ! Indices of R_{0i0j}
           ii = ivals(i); jj = jvals(i)
           !
           ! Time derivative: \pd_0 K_{m,n}
           ! NOTE we're looping over (ivals,jvals) above, which is ij ordering of indices
           !       this coincides with the ordering of first dim of kij() = xx,xy,xz,yy,yz,zz,
           !       so we can just index kij() array with i in the below
           !
           select case(nord_dt)
           ! indices tip1,...tip4 are set above depending on time-ordering of last dimension of Kij
           !   this depends on whether we are calling this rotuine from the RT or from ricci.f90
           case(2)
               ! Use second order forward difference
               dtKij = deriv1_forward2nd(kij(i,ipos,jpos,kpos,ti),kij(i,ipos,jpos,kpos,tip1), &
               & kij(i,ipos,jpos,kpos,tip2),dtf)
           case(3)
               ! Use third order forward difference
               dtKij = deriv1_forward3rd(kij(i,ipos,jpos,kpos,ti),kij(i,ipos,jpos,kpos,tip1), &
               & kij(i,ipos,jpos,kpos,tip2),kij(i,ipos,jpos,kpos,tip3),dtf)
           case(4)
               ! Use fourth order forward difference
               dtKij = deriv1_forward4th(kij(i,ipos,jpos,kpos,ti),kij(i,ipos,jpos,kpos,tip1),&
               & kij(i,ipos,jpos,kpos,tip2),kij(i,ipos,jpos,kpos,tip3),kij(i,ipos,jpos,kpos,tip4),dtf)
           case default
               call print_error("Only nord_dt = 2,3,4 implemented. Re-set and re-run.",2,loc)
           end select

           ! 2x sum terms, initialise then loop
           dkalp_gamsum = 0.d0; KKsum = 0.d0
           do k=1,3
               Gamma = gamijk(k,ii,jj,ipos,jpos,kpos)
               !call get_christoffel(ipos,jpos,kpos,gij,nx,dx,k,ii,jj,Gamma)
               dkalp_gamsum = dkalp_gamsum + alpijk * dialp(k) * Gamma
               do l=1,3 ! we have this extra loop to avoid having to calculate kud and pass it in
                   KKsum = KKsum + alpijk**2 * kdown(k,jj) * gup(k,l) * kdown(l,ii)
               enddo
           enddo

           R0i0j(ii,jj) = alpijk * dtKij + alpijk * didjalp(ii,jj) - dkalp_gamsum + KKsum
       enddo
       !
       ! Apply symmetries for remaining 3 components
       !
       R0i0j(2,1) = R0i0j(1,2)
       R0i0j(3,1) = R0i0j(1,3)
       R0i0j(3,2) = R0i0j(2,3)

   end subroutine get_R0i0j



    !
    ! Return R_{0ijk} at a particular position and indices i,j,k = idx1,idx2,idx3
    !
    ! For this one we take indices i,j,k as input and only output a single value because
    !    for the raytracer we need 9/27 components so it's not worth calculating the full
    !     (i,j,k) tensor here
    !
    ! ** checked Mar 31st 2021 - HMac OK
    !
    subroutine get_R0ijk(nx,ipos,jpos,kpos,dx,kij,alp,gamijk,kdown,idx1,idx2,idx3,R0ijk)
        integer, intent(in) :: nx,ipos,jpos,kpos  ! spatial position
        integer, intent(in) :: idx1,idx2,idx3     ! i,j,k indices of output Riemann
        real(c_double), intent(in) :: kij(6,nx,nx,nx)
        real(c_double), intent(in) :: gamijk(3,3,3,nx,nx,nx)
        real(c_double), intent(in) :: alp(nx,nx,nx),dx,kdown(3,3)

        real(c_double), intent(out) :: R0ijk

        ! Notation of indices in these terms is consistent with the eqns in user guide
        !      and, e.g., 'gamsum1' is the FIRST Gamma sum appearing in the expr.
        real(c_double) :: alpijk,gamsum1,gamsum2,gam1,gam2,djKik,dkKij
        real(c_double), dimension(3,3) :: kdown_ip1,kdown_im1,kdown_jp1,&
        & kdown_jm1,kdown_kp1,kdown_km1
        real(c_double), dimension(3,3) :: kdown_ip2,kdown_im2,kdown_jp2,kdown_jm2,&
        & kdown_kp2,kdown_km2

        integer :: ip1,im1,jp1,jm1,kp1,km1,ip2,im2,jp2,jm2,kp2,km2
        integer :: l
        character(len=clen) :: loc
        loc = " get_R0ijk"

        !
        ! Apply periodic boundary conditions
        call apply_periodic(ipos,ip1,im1,nx)
        call apply_periodic(jpos,jp1,jm1,nx)
        call apply_periodic(kpos,kp1,km1,nx)
        call apply_periodic_fourth(ipos,ip2,im2,nx)
        call apply_periodic_fourth(jpos,jp2,jm2,nx)
        call apply_periodic_fourth(kpos,kp2,km2,nx)

        !
        ! Get lapse, K_ij at current position + spatial stencil for Kij
        alpijk = alp(ipos,jpos,kpos)
        !call get_metric_at_pos(ipos,jpos,kpos,nx,kij,kdown)
        ! NOTE we could pass in kdown_atstencil here but this is just indexing...
        !    so it's probably the same as indexing kdown_atstencil here?
        call get_metric_at_pos(ip1,jpos,kpos,nx,kij,kdown_ip1)
        call get_metric_at_pos(im1,jpos,kpos,nx,kij,kdown_im1)
        call get_metric_at_pos(ipos,jp1,kpos,nx,kij,kdown_jp1)
        call get_metric_at_pos(ipos,jm1,kpos,nx,kij,kdown_jm1)
        call get_metric_at_pos(ipos,jpos,kp1,nx,kij,kdown_kp1)
        call get_metric_at_pos(ipos,jpos,km1,nx,kij,kdown_km1)

        call get_metric_at_pos(ip2,jpos,kpos,nx,kij,kdown_ip2)
        call get_metric_at_pos(im2,jpos,kpos,nx,kij,kdown_im2)
        call get_metric_at_pos(ipos,jp2,kpos,nx,kij,kdown_jp2)
        call get_metric_at_pos(ipos,jm2,kpos,nx,kij,kdown_jm2)
        call get_metric_at_pos(ipos,jpos,kp2,nx,kij,kdown_kp2)
        call get_metric_at_pos(ipos,jpos,km2,nx,kij,kdown_km2)

        !
        ! Calculate derivatives of K_ij
        !     which derivs we take are dependent on the value of indices passed in
        !
        ! \pd_j K_ik = \pd_idx2 K_{idx1,idx3}
        select case(idx2)
        case(1)
            ! x-derivative of K_{idx1,idx3}
            djKik = deriv1fourth(kdown_ip1(idx1,idx3),kdown_ip2(idx1,idx3),&
            & kdown_im1(idx1,idx3),kdown_im2(idx1,idx3),dx)
        case(2)
            ! y-derivative of K_{idx1,idx3}
            djKik = deriv1fourth(kdown_jp1(idx1,idx3),kdown_jp2(idx1,idx3),&
            & kdown_jm1(idx1,idx3),kdown_jm2(idx1,idx3),dx)
        case(3)
            ! z-derivative of K_{idx1,idx3}
            djKik = deriv1fourth(kdown_kp1(idx1,idx3),kdown_kp2(idx1,idx3),&
            & kdown_km1(idx1,idx3),kdown_km2(idx1,idx3),dx)
        case default
            ! weird value sent in
            djKik = -999
            call print_error("idx2 not 1,2 or 3. Check yo' self.",2,loc)
        end select

        ! \pd_k K_{ij} = \pd_idx3 K_{idx1,idx2}
        select case(idx3)
        case(1)
            ! x-derivative of K_{idx1,idx2}
            dkKij = deriv1fourth(kdown_ip1(idx1,idx2),kdown_ip2(idx1,idx2),&
            & kdown_im1(idx1,idx2),kdown_im2(idx1,idx2),dx)
        case(2)
            ! y-derivative of K_{idx1,idx2}
            dkKij = deriv1fourth(kdown_jp1(idx1,idx2),kdown_jp2(idx1,idx2),&
            & kdown_jm1(idx1,idx2),kdown_jm2(idx1,idx2),dx)
        case(3)
            ! z-derivative of K_{idx1,idx2}
            dkKij = deriv1fourth(kdown_kp1(idx1,idx2),kdown_kp2(idx1,idx2),&
            & kdown_km1(idx1,idx2),kdown_km2(idx1,idx2),dx)
        case default
            ! weird value sent in
            dkKij = -999
            call print_error("idx2 not 1,2 or 3. Check yo' self.",2,loc)
        end select

        gamsum1 = 0.d0; gamsum2 = 0.d0
        do l=1,3
            ! alp K_lj \Gam^l_{ik} = alp K_{l,idx2} \Gam^l_{idx1,idx3}
            gam1 = gamijk(l,idx1,idx3,ipos,jpos,kpos)
            !call get_christoffel(ipos,jpos,kpos,gij,nx,dx,l,idx1,idx3,gam1)
            gamsum1 = gamsum1 + alpijk * kdown(l,idx2) * gam1

            ! alp K_lk \Gam^l_{ij} = alp K_{l,idx3} \Gam^l_{idx1,idx2}
            gam2 = gamijk(l,idx1,idx2,ipos,jpos,kpos)
            !call get_christoffel(ipos,jpos,kpos,gij,nx,dx,l,idx1,idx2,gam2)
            gamsum2 = gamsum2 + alpijk * kdown(l,idx3) * gam2
        enddo

        R0ijk = alpijk * djKik - alpijk * dkKij + gamsum1 - gamsum2

    end subroutine get_R0ijk



    !
    ! Return R_{ijkl} at a particular position and indices i,j,k,l = idx1,idx2,idx3,idx4
    !
    ! For this one we take indices i,j,k as input and only output a single value because
    !    for the raytracer we need 6/(3x3x3x3) components so it's not worth calculating the full
    !     (i,j,k,l) tensor here
    !
    ! NOTE: 23/07/2021: WAY slower to call get_christoffel in here (and all thru RT)
    !            trust past Hay... don't do it
    !
    subroutine get_Rijkl(nx,dx,ipos,jpos,kpos,gamijk,gdown,kdown,idx1,idx2,idx3,idx4,Rijkl)
        integer, intent(in) :: nx,ipos,jpos,kpos       ! spatial position
        integer, intent(in) :: idx1,idx2,idx3,idx4     ! i,j,k,l indices of output Riemann
        real(c_double), intent(in) :: gamijk(3,3,3,nx,nx,nx)
        real(c_double), intent(in) :: dx,gdown(3,3),kdown(3,3)

        real(c_double), intent(out) :: Rijkl

        ! Index Notation here is consistent with form of eqns in user guide
        !    and, e.g., GamGam1 is the first double Gamma sum term that appears in the expr.
        real(c_double) :: dkGam_mjl,dlGam_mjk,GamGam1,GamGam2
        ! Christoffels at stencil for derivative terms above
        ! \Gam^m_{jl} and \Gam^m_{jk} at stencil, depending which derivative we take
        real(c_double) :: Gammjl_p1,Gammjl_p2,Gammjl_m1,Gammjl_m2
        real(c_double) :: Gammjk_p1,Gammjk_p2,Gammjk_m1,Gammjk_m2
        real(c_double) :: gam_mnk,gam_njl,gam_mnl,gam_njk

        integer :: ip1,im1,jp1,jm1,kp1,km1,ip2,im2,jp2,jm2,kp2,km2 ! for periodic boundaries
        integer :: m,n
        character(len=clen) :: loc
        loc = " get_Rijkl"

        !
        ! Apply periodic boundary conditions
        call apply_periodic(ipos,ip1,im1,nx)
        call apply_periodic(jpos,jp1,jm1,nx)
        call apply_periodic(kpos,kp1,km1,nx)
        call apply_periodic_fourth(ipos,ip2,im2,nx)
        call apply_periodic_fourth(jpos,jp2,jm2,nx)
        call apply_periodic_fourth(kpos,kp2,km2,nx)

        !
        ! Zero-sum terms: Kij's
        Rijkl = kdown(idx1,idx3) * kdown(idx2,idx4) - kdown(idx1,idx4) * kdown(idx2,idx3)
        do m=1,3
            !
            ! First, get derivative terms of Christoffels
            !   this will depend on which indices are passed in
            !
            ! \pd_k \Gam^m_{jl} = \pd_idx3 \Gam^dum_{idx2,idx4}
            select case(idx3)
            case(1)
                ! x-derivative of \Gam^dum_{idx2,idx4}
                !
                ! 1. get christoffel at x-stencil
                Gammjl_p1 = gamijk(m,idx2,idx4,ip1,jpos,kpos)
                Gammjl_p2 = gamijk(m,idx2,idx4,ip2,jpos,kpos)
                Gammjl_m1 = gamijk(m,idx2,idx4,im1,jpos,kpos)
                Gammjl_m2 = gamijk(m,idx2,idx4,im2,jpos,kpos)
                !call get_christoffel(ip1,jpos,kpos,gij,nx,dx,m,idx2,idx4,Gammjl_p1)
                !call get_christoffel(ip2,jpos,kpos,gij,nx,dx,m,idx2,idx4,Gammjl_p2)
                !call get_christoffel(im1,jpos,kpos,gij,nx,dx,m,idx2,idx4,Gammjl_m1)
                !call get_christoffel(im2,jpos,kpos,gij,nx,dx,m,idx2,idx4,Gammjl_m2)
            case(2)
                ! y-derivative of \Gam^dum_{idx2,idx4}
                !
                ! 1. get christoffel at y-stencil
                Gammjl_p1 = gamijk(m,idx2,idx4,ipos,jp1,kpos)
                Gammjl_p2 = gamijk(m,idx2,idx4,ipos,jp2,kpos)
                Gammjl_m1 = gamijk(m,idx2,idx4,ipos,jm1,kpos)
                Gammjl_m2 = gamijk(m,idx2,idx4,ipos,jm2,kpos)
                !call get_christoffel(ipos,jp1,kpos,gij,nx,dx,m,idx2,idx4,Gammjl_p1)
                !call get_christoffel(ipos,jp2,kpos,gij,nx,dx,m,idx2,idx4,Gammjl_p2)
                !call get_christoffel(ipos,jm1,kpos,gij,nx,dx,m,idx2,idx4,Gammjl_m1)
                !call get_christoffel(ipos,jm2,kpos,gij,nx,dx,m,idx2,idx4,Gammjl_m2)
            case(3)
                ! z-derivative of \Gam^dum_{idx2,idx4}
                !
                ! 1. get christoffel at z-stencil
                Gammjl_p1 = gamijk(m,idx2,idx4,ipos,jpos,kp1)
                Gammjl_p2 = gamijk(m,idx2,idx4,ipos,jpos,kp2)
                Gammjl_m1 = gamijk(m,idx2,idx4,ipos,jpos,km1)
                Gammjl_m2 = gamijk(m,idx2,idx4,ipos,jpos,km2)
                !call get_christoffel(ipos,jpos,kp1,gij,nx,dx,m,idx2,idx4,Gammjl_p1)
                !call get_christoffel(ipos,jpos,kp2,gij,nx,dx,m,idx2,idx4,Gammjl_p2)
                !call get_christoffel(ipos,jpos,km1,gij,nx,dx,m,idx2,idx4,Gammjl_m1)
                !call get_christoffel(ipos,jpos,km2,gij,nx,dx,m,idx2,idx4,Gammjl_m2)
            end select
            ! 2. take derivative
            dkGam_mjl = deriv1fourth(Gammjl_p1,Gammjl_p2,Gammjl_m1,Gammjl_m2,dx)

            ! \pd_l \Gam^m_{jk} = \pd_idx4 \Gam^dum_{idx2,idx3}
            select case(idx4)
            case(1)
                ! x-derivative of \Gam^dum_{idx2,idx3}
                !
                ! 1. get christoffel at x-stencil
                Gammjk_p1 = gamijk(m,idx2,idx3,ip1,jpos,kpos)
                Gammjk_p2 = gamijk(m,idx2,idx3,ip2,jpos,kpos)
                Gammjk_m1 = gamijk(m,idx2,idx3,im1,jpos,kpos)
                Gammjk_m2 = gamijk(m,idx2,idx3,im2,jpos,kpos)
                !call get_christoffel(ip1,jpos,kpos,gij,nx,dx,m,idx2,idx3,Gammjk_p1)
                !call get_christoffel(ip2,jpos,kpos,gij,nx,dx,m,idx2,idx3,Gammjk_p2)
                !call get_christoffel(im1,jpos,kpos,gij,nx,dx,m,idx2,idx3,Gammjk_m1)
                !call get_christoffel(im2,jpos,kpos,gij,nx,dx,m,idx2,idx3,Gammjk_m2)
            case(2)
                ! y-derivative of \Gam^dum_{idx2,idx3}
                !
                ! 1. get christoffel at y-stencil
                Gammjk_p1 = gamijk(m,idx2,idx3,ipos,jp1,kpos)
                Gammjk_p2 = gamijk(m,idx2,idx3,ipos,jp2,kpos)
                Gammjk_m1 = gamijk(m,idx2,idx3,ipos,jm1,kpos)
                Gammjk_m2 = gamijk(m,idx2,idx3,ipos,jm2,kpos)
                !call get_christoffel(ipos,jp1,kpos,gij,nx,dx,m,idx2,idx3,Gammjk_p1)
                !call get_christoffel(ipos,jp2,kpos,gij,nx,dx,m,idx2,idx3,Gammjk_p2)
                !call get_christoffel(ipos,jm1,kpos,gij,nx,dx,m,idx2,idx3,Gammjk_m1)
                !call get_christoffel(ipos,jm2,kpos,gij,nx,dx,m,idx2,idx3,Gammjk_m2)
            case(3)
                ! z-derivative of \Gam^dum_{idx2,idx3}
                !
                ! 1. get christoffel at z-stencil
                Gammjk_p1 = gamijk(m,idx2,idx3,ipos,jpos,kp1)
                Gammjk_p2 = gamijk(m,idx2,idx3,ipos,jpos,kp2)
                Gammjk_m1 = gamijk(m,idx2,idx3,ipos,jpos,km1)
                Gammjk_m2 = gamijk(m,idx2,idx3,ipos,jpos,km2)
                !call get_christoffel(ipos,jpos,kp1,gij,nx,dx,m,idx2,idx3,Gammjk_p1)
                !call get_christoffel(ipos,jpos,kp2,gij,nx,dx,m,idx2,idx3,Gammjk_p2)
                !call get_christoffel(ipos,jpos,km1,gij,nx,dx,m,idx2,idx3,Gammjk_m1)
                !call get_christoffel(ipos,jpos,km2,gij,nx,dx,m,idx2,idx3,Gammjk_m2)
            end select
            ! 2. take derivative
            dlGam_mjk = deriv1fourth(Gammjk_p1,Gammjk_p2,Gammjk_m1,Gammjk_m2,dx)

            !
            ! Double gamma terms -- double sum
            GamGam1 = 0.d0; GamGam2 = 0.d0
            do n=1,3
                ! \Gam^m_{nk} \Gam^n_{jl} = \Gam^m_{n,idx3} \Gam^n_{idx2,idx4}
                gam_mnk = gamijk(m,n,idx3,ipos,jpos,kpos)
                !call get_christoffel(ipos,jpos,kpos,gij,nx,dx,m,n,idx3,gam_mnk)
                gam_njl = gamijk(n,idx2,idx4,ipos,jpos,kpos)
                !call get_christoffel(ipos,jpos,kpos,gij,nx,dx,n,idx2,idx4,gam_njl)
                GamGam1 = GamGam1 + gam_mnk * gam_njl

                ! \Gam^m_{nl} \Gam^n_{jk} = \Gam^m_{n,idx4} \Gam^n_{idx2,idx3}
                gam_mnl = gamijk(m,n,idx4,ipos,jpos,kpos)
                !call get_christoffel(ipos,jpos,kpos,gij,nx,dx,m,n,idx4,gam_mnl)
                gam_njk = gamijk(n,idx2,idx3,ipos,jpos,kpos)
                !call get_christoffel(ipos,jpos,kpos,gij,nx,dx,n,idx2,idx3,gam_njk)
                GamGam2 = GamGam2 + gam_mnl * gam_njk
            enddo

            ! sum over m index included in here
            Rijkl = Rijkl + gdown(idx1,m) * (dkGam_mjl - dlGam_mjk + GamGam1 &
                & - GamGam2)

        enddo

    end subroutine get_Rijkl




end module riemann
