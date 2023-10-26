module weyl
    !
    ! A module to contain routines related to calculating
    !   anything related to the Weyl tensor in regular Mescaline
    !    -- i.e., to be called from Ricci (not RT...)
    !
    use options, only:c_double,clen
    use riemann, only:get_RAB
    use prints, only:print_error
    implicit none

contains

    !
    ! A subroutine to return the ELECTRIC part of the Weyl tensor
    ! -- see written notes:
    !
    ! E_{alp,beta} = 4 C_{AB} F^[A]_alp F^[B]_beta
    !    where F^A_alp \equiv u^mu h^nu_alp and [] is antisymmetric part
    !
    !    -- checked 31/08/2021 -- OK HMAC
    subroutine get_electric_weyl(nx,nt,ti,ipos,jpos,kpos,dx,dt,gij,kij,alp,gamijk,gdown,&
        & gup,kdown,fourR,Rdd,umuU,umud,Edd,Esq)
        integer, intent(in) :: nx,nt,ti,ipos,jpos,kpos
        real(c_double), intent(in) :: dx,dt,fourR,Rdd(4,4)
        real(c_double), intent(in) :: umuU(4),umud(4) ! u^mu and u_mu (both calculated in get_expansion_shear_vort)
        real(c_double), dimension(3,3), intent(in) :: gdown,gup,kdown
        real(c_double), dimension(6,nx,nx,nx,nt), intent(in) :: kij
        real(c_double), dimension(6,nx,nx,nx),    intent(in) :: gij
        real(c_double), intent(in) :: alp(nx,nx,nx),gamijk(3,3,3,nx,nx,nx)

        real(c_double), intent(out) :: Edd(4,4)      ! Electric Weyl
        real(c_double), intent(out), optional :: Esq ! + scalar Esq = EUU Edd / 2

        real(c_double) :: CAB(6,6),FF(6,4),FB(6,4),Fas(6,4) ! the Weyl tensor, F^Af_alp, F^Ab_alp, and F^[A]_alp
        real(c_double) :: deltamunu,hud(4,4),gUUmunu(4,4) ! projection tensor h^alp_beta
        real(c_double) :: trE

        integer :: mu,nu,alpha,bet,A,B
        logical :: checktrE

        ! 1. Get the Weyl tensor, CAB(6,6)
        call get_CAB(nx,nt,ti,ipos,jpos,kpos,dx,dt,gij,kij,alp,gamijk,gdown,&
            & gup,kdown,fourR,Rdd,CAB)

        !
        ! 2. Set up F^[A]_alp = 0.5 * ( F^Af_alp - F^Ab_alp ) = 0.5 * (FF - FB)
        !
        !     a) set up hud = h^alp_beta = g^alp_beta + u^alp u_beta = delta^alp_beta + u^alp u_beta
        do mu=1,4
            do nu=1,4
                if (mu==nu) then
                    ! delta^mu_nu = 1
                    deltamunu = 1.d0
                else
                    ! delta^mu_nu = 0
                    deltamunu = 0.d0
                endif
                hud(mu,nu) = deltamunu + umuU(mu) * umud(nu)
            enddo
        enddo
        ! Now F^Af_alp = u^mu h^nu_alp for Af = (mu,nu) = (01,02,03,12,13,23)
        !  just do this by hand, seems easiest...

        ! Af=1 --> mu,nu = 0,1 (index 1,2)
        FF(1,:) = umuU(1) * hud(2,:)
        ! Af=2 --> mu,nu = 0,2 (index 1,3)
        FF(2,:) = umuU(1) * hud(3,:)
        ! Af=3 --> mu,nu = 0,3 (index 1,4)
        FF(3,:) = umuU(1) * hud(4,:)
        ! Af=4 --> mu,nu = 1,2 (index 2,3)
        FF(4,:) = umuU(2) * hud(3,:)
        ! Af=5 --> mu,nu = 1,3 (index 2,4)
        FF(5,:) = umuU(2) * hud(4,:)
        ! Af=6 --> mu,nu = 2,3 (index 3,4)
        FF(6,:) = umuU(3) * hud(4,:)

        ! Now F^Ab_alp = u^mu h^nu_alp for Ab = (mu,nu) = (10,20,30,21,31,32)

        ! Ab=1 --> mu,nu = 1,0 (index 2,1)
        FB(1,:) = umuU(2) * hud(1,:)
        ! Ab=2 --> mu,nu = 2,0 (index 3,1)
        FB(2,:) = umuU(3) * hud(1,:)
        ! Ab=3 --> mu,nu = 3,0 (index 4,1)
        FB(3,:) = umuU(4) * hud(1,:)
        ! Ab=4 --> mu,nu = 2,1 (index 3,2)
        FB(4,:) = umuU(3) * hud(2,:)
        ! Ab=5 --> mu,nu = 3,1 (index 4,2)
        FB(5,:) = umuU(4) * hud(2,:)
        ! Ab=6 --> mu,nu = 3,2 (index 4,3)
        FB(6,:) = umuU(4) * hud(3,:)

        ! Now we can construct the antisymmetric part, F^[A]_alp = F^Af_alp - F^Ab_alp
        !   (we neglect the 1/2 factor, because there is a 4* outside sum in E_{alp,beta}
        !        that cancels when we consider these.)
        Fas = FF - FB

        ! 3. Compute E_{alp,bet} by looping over components and summing (A,B) over 6,6 loop
        Edd = 0.d0
        do alpha=1,4
            do bet=1,4

                ! sum loop over (A,B)
                do A=1,6
                    do B=1,6
                        Edd(alpha,bet) = Edd(alpha,bet) + CAB(A,B) * Fas(A,alpha) * Fas(B,bet)
                        !if (abs(CAB(A,B))>1.e-15) then
                        !    print*, 'A,B, CAB(A,B) = ',A,B,CAB(A,B)
                        !endif
                    enddo
                enddo

            enddo
        enddo

        if (present(Esq)) then
            ! get Esq scalar as well
            ! Scalar Esq = E^2 = E^{ab} E_{ab} = g^{ac} g^{bd} E_{cd} E_{ab} = 4x sum, if A,B<5
            !
            ! 1. build g^{mu,nu} for Esq sum
            gUUmunu      = 0.d0
            gUUmunu(1,1) = -1.d0 / alp(ipos,jpos,kpos)**2 ! g^{00}
            gUUmunu(2:4,2:4) = gup
            Esq = 0.d0; trE = 0.d0
            do alpha=1,4
                do bet=1,4
                    do A=1,4
                        do B=1,4
                            Esq = Esq + 0.5d0 * gUUmunu(alpha,A) * gUUmunu(bet,B) * Edd(A,B) * Edd(alpha,bet)
                        enddo
                    enddo
                    !
                    ! calculate trace of E_mu,nu. Should be zero; so a good check of errors
                    trE = trE + gUUmunu(alpha,bet) * Edd(alpha,bet)
                enddo
            enddo
        endif
        checktrE = .False.
        if (checktrE) then
           ! check trE is small
           if (abs(trE)>1.e-14) print*, 'Large trace E_{mu,nu}: trE = ',trE,' E_00, E_0x, E_xx, E_xy = ',&
                & Edd(1,1), Edd(1,2), Edd(2,2), Edd(2,3)
        endif


    end subroutine get_electric_weyl





    !
    ! A subroutine to return the components of the Weyl tensor
    !     C_{AB} where A,B\in {01,02,03,12,13,23}
    !
    !    -- checked 31/08/2021 -- OK HMAC
    subroutine get_CAB(nx,nt,ti,ipos,jpos,kpos,dx,dt,gij,kij,alp,gamijk,gdown,&
        & gup,kdown,fourR,Rdd,CAB)
        integer, intent(in) :: nx,nt,ti,ipos,jpos,kpos
        real(c_double), intent(in) :: dx,dt,fourR,Rdd(4,4)
        real(c_double), dimension(3,3), intent(in) :: gdown,gup,kdown
        real(c_double), dimension(6,nx,nx,nx,nt), intent(in) :: kij
        real(c_double), dimension(6,nx,nx,nx),    intent(in) :: gij
        real(c_double), intent(in) :: alp(nx,nx,nx),gamijk(3,3,3,nx,nx,nx)

        real(c_double), intent(out) :: CAB(6,6) ! the Weyl tensor

        real(c_double) :: RAB(6,6),gmunu(4,4),ricci_part
        integer :: i,j,mus(6),nus(6),mu,nu,lam,rho
        logical :: check

        ! 1. First step; get the Riemann tensor R_{AB}

        !
        ! TODO: you will need to somehow get around the time ordering of kij array
        !    --> in get_R0i0j we need dt Kij, and this is currently coded for the
        !         ordering of time shots for the raytracer, which is opposite to regular mescaline
        !    --> maybe add a flag that can be passed in, or change one or the other so that they're the same?
        call get_RAB(nx,nt,ti,ipos,jpos,kpos,dx,dt,gij,kij,alp,gamijk,gdown,gup,kdown,RAB)

        ! 2. Build the 4D tensor gmunu so we can easily pass it into get_ricci_part_of_weyl for indexing
        gmunu      = 0.d0
        ! a) g_00 = - alp^2
        gmunu(1,1) = - alp(ipos,jpos,kpos)**2
        ! b) g_ij = gdown (gamma_ij)
        gmunu(2:4,2:4) = gdown

        ! 3. Set up the indices so we can loop over Weyl components
        !   C_{AB} = C_{mu,nu,lambda,rho}
        !
        !   the indices A,B \in [01,02,03,12,13,23]

        !   this means (mu,lambda) \in [0,0,0,1,1,2]
        !        and      (nu,rho) \in [1,2,3,2,3,3]
        mus = (/ 0, 0, 0, 1, 1, 2 /) ! and lambdas
        nus = (/ 1, 2, 3, 2, 3, 3 /) ! and rhos

        ! Loop over A,B
        do i=1,6
            do j=1,6
                ! Set up mu,nu,lam,rho indices for this value of A,B
                ! Actual indices need to +1 for F90
                mu  = mus(i)+1;  nu = nus(i)+1    ! value of A = (mu,nu)
                lam = mus(j)+1; rho = nus(j)+1  ! value of B = (lam,rho)

                ! Get the remaining part of the Weyl tensor (non-Riemann tensor part = Ricci part)
                call get_ricci_part_of_weyl(mu,nu,lam,rho,gmunu,Rdd,fourR,ricci_part)

                CAB(i,j) = RAB(i,j) + ricci_part
            enddo
        enddo

        ! Do a sanity check that contracting either A or B gives zero
        check = .True.
        if (check) call check_weyl_contraction(CAB,alp(ipos,jpos,kpos),gup)

    end subroutine get_CAB


    !
    ! A subroutine to return the "Ricci part" of the Weyl tensor
    !     namely, the part that is not the Riemann tensor + contains the 4-Ricci tensor & scalar
    !  --> C_AB = R_AB + ricci_part (where AB = mu,nu,lambda,rho)
    !
    !    -- checked 31/08/2021 -- OK HMAC
    subroutine get_ricci_part_of_weyl(mu,nu,lam,rho,gdd,Rdd,R,ricci_part)
        integer, intent(in) :: mu,nu,lam,rho ! the indices of the Weyl tensor component we are calculating
        real(c_double), intent(in) :: gdd(4,4),Rdd(4,4),R ! g_mu,nu, 4-Ricci tensor & scalar

        real(c_double), intent(out) :: ricci_part

        real(c_double) :: term1,term2

        ! Term with 1/2 factor containing both g_munu and R_munu
        term1 = gdd(mu,rho) * Rdd(nu,lam) + gdd(nu,lam) * Rdd(mu,rho) &
            & - gdd(mu,lam) * Rdd(nu,rho) - gdd(nu,rho) * Rdd(mu,lam)

        ! Term multiplied by the Ricci scalar containing only g_munu
        term2 = gdd(mu,lam) * gdd(nu,rho) - gdd(mu,rho) * gdd(nu,lam)

        ricci_part = 0.5d0 * term1 + R * term2 / 6.d0

    end subroutine get_ricci_part_of_weyl





    !
    ! A subroutine to check that raising + contracting some test indices
    ! of C_AB gives zero (this should always be true)
    !
    subroutine check_weyl_contraction(CAB,alp,gup)
        real(c_double), intent(in) :: CAB(6,6),alp,gup(3,3)

        real(c_double) :: gupmunu(4,4),CAB_Asum(6),CAB_Bsum(6)
        integer :: i,j,mus(6),nus(6),mu,nu,lam,rho
        character(len=clen), parameter :: loc=" check_weyl_tracefree_index1"
        ! g^{mu,nu} C_{mu,nu,lam,rho} = g^A C_{AB} = 0
        ! g^{lam,rho} C_{mu,nu,lam,rho} = g^B C_{AB} = 0

        gupmunu = 0.d0
        gupmunu(1,1) = -1.d0 / alp**2
        gupmunu(2:4,2:4) = gup

        mus = (/ 0, 0, 0, 1, 1, 2 /) ! and lambdas
        nus = (/ 1, 2, 3, 2, 3, 3 /) ! and rhos

        ! Loop over A,B
        CAB_Asum = 0.d0; CAB_Bsum = 0.d0
        do i=1,6 ! A
            do j=1,6
                ! B
                mu  = mus(i)+1;  nu = nus(i)+1    ! value of A = (mu,nu)
                lam = mus(j)+1; rho = nus(j)+1    ! value of B = (lam,rho)
                CAB_Asum(j) = CAB_Asum(j) + gupmunu(mu,nu) * CAB(i,j)
                CAB_Bsum(i) = CAB_Bsum(i) + gupmunu(lam,rho) * CAB(i,j)
            enddo
        enddo

        ! check if all components are zero
        do j=1,6
           if (abs(CAB_Asum(j))>1.e-14) then
              call print_error("g^A C_{AB} > 1.e-15",0,loc)
              print*, 'j,CAB_Asum(j) = ',j,CAB_Asum(j)
           endif
           if (abs(CAB_Bsum(j))>1.e-14) then
              call print_error("g^B C_{AB} > 1.e-15",0,loc)
              print*, 'j,CAB_Bsum(j) = ',j,CAB_Bsum(j)
           endif
           !print*, 'j,CAB_Asum(j), CAB_Bsum(j) = ',j,CAB_Asum(j),CAB_Bsum(j)

        enddo

    end subroutine check_weyl_contraction




end module weyl
