module violation
  !----------------------------------------------
  !
  ! Module to compute the level of Hamiltonian (H) and
  !     momentum (M_i) violation including relative energy scale
  !      for normalisation
  !
  !----------------------------------------------
  use periodic, only: apply_periodic,apply_periodic_fourth
  use manipulations, only:get_metric_at_pos,inv3x3,trace,calc_up_down,calc_down_up
  use derivatives, only:deriv1fourth
  use tripreports, only:filename_simp,write_avg,write_3d_data
  use options, only:dp,pi,c_double,rawHM,rawandrelHM,clen
  use prints, only:print_info
  implicit none
contains


  !
  ! subroutine to calculate Hamiltonian constraint violation IN main loop (i.e. called at every i,j,k)
  !
  subroutine calc_ham_violation(tracer,tracek,kud,rho,kudud,violation,e_scale2)
    real(c_double), intent(in) :: tracer, tracek, rho
    real(c_double), intent(in), dimension(3,3) :: kud
    real(c_double), intent(out) :: kudud,e_scale2, violation ! relative Ham violation (H/[H])
    integer :: l,m

    kudud = 0._dp
    do l=1,3
       do m=1,3
          kudud = kudud + kud(l,m) * kud(m,l)
       enddo
    enddo
    !
    ! calculate "energy scale" [H]
    e_scale2  = tracer**2 + tracek**4 + kudud**2 + (16._dp * pi*rho)**2
    !
    ! Hamiltonian constraint violation
    violation = tracer - kudud + tracek**2 - 16._dp * pi*rho

  end subroutine calc_ham_violation




  !
  ! Calculate the constraint violation (and write to file)
  ! relative to energy scale:
  !      M_i = D_i K^j_i - D_j trK - S_i --> and M = sqrt{M^i M_i}
  !             (D_i is cov deriv associated with \gamma_ij)
  !        H = 3R - K_ij K^ij - K^2 + 16pirho
  !
  ! NOTE: this routine is to be called from WITHIN a spatial loop. this is called from the main ricci loop ONLY
  !
  ! -- HJM 29/09/2023: checked all Mom terms + re-derived term1 - all OK!
  !                    checked Christoffel symmetries satisfied exactly in mesc test - OK
  !
  subroutine constraint_violation_inloop(i,j,k,nx,dx,gij_atstencil,kij_atstencil,Chrsijk,&
       tracek,rho,vel1,vel2,vel3,lorentz,tracer,rhoAij,Add,AUU,magMom2,Momescale2,Hamraw,Hamescale2)
    integer, intent(in) :: i,j,k,nx
    real(c_double), intent(in) :: dx,Chrsijk(3,3,3,nx,nx,nx)  ! need derivatives of this to get 3R
    real(c_double), intent(in), dimension(nx,nx,nx) :: tracek      ! need derivatives of this for M_i
    real(c_double), intent(in) :: rho,vel1,vel2,vel3,lorentz,tracer
    !real(c_double), intent(in), dimension(6,nx,nx,nx) :: gij
    real(c_double), intent(in), dimension(13,3,3) :: gij_atstencil,kij_atstencil

    real(c_double), intent(out) :: Add(3,3),AUU(3,3),rhoAij ! A_ij and rho_{Aij} -- amplitude of tensor inhomog. (Clough+2018)
    real(c_double), intent(out) :: magMom2,Momescale2      ! Output |M|^2 to avoid sqrt here
    real(c_double), intent(out) :: Hamraw,Hamescale2       ! e_scales^2 as well for same reason

    real(c_double) :: term1a,term1b,term1c
    real(c_double) :: term1mag,term2mag,term3mag
    real(c_double) :: gamma1b_i,gamma1b_ii,gamma1b_iii,gamma1c_i,gamma1c_ii,gamma1c_iii
    real(c_double) :: trKijk,trKip1,trKim1,trKjp1,trKjm1,trKkp1,trKkm1,detgijk,dum
    real(c_double) :: trKip2,trKim2,trKjp2,trKjm2,trKkp2,trKkm2,kudud,kup(3,3)
    real(c_double) :: rho_n

    real(c_double), dimension(3)   :: Mom,term1,term2,term3
    real(c_double), dimension(3,3) :: kud,kudip1,kudim1,kudjp1,kudjm1,kudkp1,kudkm1
    real(c_double), dimension(3,3) :: kudip2,kudim2,kudjp2,kudjm2,kudkp2,kudkm2
    real(c_double), dimension(3,3) :: gdown,gdownip1,gdownim1,gdownjp1,gdownjm1,gdownkp1,gdownkm1
    real(c_double), dimension(3,3) :: gdownip2,gdownim2,gdownjp2,gdownjm2,gdownkp2,gdownkm2
    real(c_double), dimension(3,3) :: gup,gupip1,gupim1,gupjp1,gupjm1,gupkp1,gupkm1
    real(c_double), dimension(3,3) :: gupip2,gupim2,gupkp2,gupkm2,gupjp2,gupjm2
    real(c_double), dimension(3,3) :: kdown,kdownip1,kdownim1,kdownjp1,kdownjm1,kdownkp1,kdownkm1
    real(c_double), dimension(3,3) :: kdownip2,kdownim2,kdownjp2,kdownjm2,kdownkp2,kdownkm2

    integer :: ip1,jp1,kp1,im1,jm1,km1,l,m,n
    integer :: ip2,im2,jp2,jm2,kp2,km2

    !
    ! Get K_ij at all positions
    !
    !      --> gij_atstencil first dim is (pos,im2,im1,ip1,ip2,jm2,jm1,jp1,jp2,km2,km1,kp1,kp2)
    !
    kdown = kij_atstencil(1,:,:)
    kdownim2 = kij_atstencil(2,:,:)
    kdownip2 = kij_atstencil(5,:,:)
    kdownjm2 = kij_atstencil(6,:,:)
    kdownjp2 = kij_atstencil(9,:,:)
    kdownkm2 = kij_atstencil(10,:,:)
    kdownkp2 = kij_atstencil(13,:,:)

    kdownim1 = kij_atstencil(3,:,:)
    kdownip1 = kij_atstencil(4,:,:)
    kdownjm1 = kij_atstencil(7,:,:)
    kdownjp1 = kij_atstencil(8,:,:)
    kdownkm1 = kij_atstencil(11,:,:)
    kdownkp1 = kij_atstencil(12,:,:)

    !
    ! Get g_ij and g^ij at all positions
    !
    !       --> gij_atstencil first dim is (pos,im2,im1,ip1,ip2,jm2,jm1,jp1,jp2,km2,km1,kp1,kp2)
    !
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

    call inv3x3(gdown,gup,detgijk)
    call inv3x3(gdownip1,gupip1,dum)
    call inv3x3(gdownim1,gupim1,dum)
    call inv3x3(gdownjp1,gupjp1,dum)
    call inv3x3(gdownjm1,gupjm1,dum)
    call inv3x3(gdownkp1,gupkp1,dum)
    call inv3x3(gdownkm1,gupkm1,dum)

    call inv3x3(gdownip2,gupip2,dum)
    call inv3x3(gdownim2,gupim2,dum)
    call inv3x3(gdownjp2,gupjp2,dum)
    call inv3x3(gdownjm2,gupjm2,dum)
    call inv3x3(gdownkp2,gupkp2,dum)
    call inv3x3(gdownkm2,gupkm2,dum)

    ! get K^i_j at all positions
    call calc_up_down(kdown,gup,kud)
    call calc_up_down(kdownip1,gupip1,kudip1)
    call calc_up_down(kdownim1,gupim1,kudim1)
    call calc_up_down(kdownjp1,gupjp1,kudjp1)
    call calc_up_down(kdownjm1,gupjm1,kudjm1)
    call calc_up_down(kdownkp1,gupkp1,kudkp1)
    call calc_up_down(kdownkm1,gupkm1,kudkm1)

    call calc_up_down(kdownip2,gupip2,kudip2)
    call calc_up_down(kdownim2,gupim2,kudim2)
    call calc_up_down(kdownjp2,gupjp2,kudjp2)
    call calc_up_down(kdownjm2,gupjm2,kudjm2)
    call calc_up_down(kdownkp2,gupkp2,kudkp2)
    call calc_up_down(kdownkm2,gupkm2,kudkm2)

    !
    ! Get traceK at all positions
    !
    call apply_periodic(i,ip1,im1,nx)
    call apply_periodic(j,jp1,jm1,nx)
    call apply_periodic(k,kp1,km1,nx)
    call apply_periodic_fourth(i,ip2,im2,nx)
    call apply_periodic_fourth(j,jp2,jm2,nx)
    call apply_periodic_fourth(k,kp2,km2,nx)

    trKip1 = tracek(ip1,j,k)
    trKim1 = tracek(im1,j,k)
    trKjp1 = tracek(i,jp1,k)
    trKjm1 = tracek(i,jm1,k)
    trKkp1 = tracek(i,j,kp1)
    trKkm1 = tracek(i,j,km1)

    trKip2 = tracek(ip2,j,k)
    trKim2 = tracek(im2,j,k)
    trKjp2 = tracek(i,jp2,k)
    trKjm2 = tracek(i,jm2,k)
    trKkp2 = tracek(i,j,kp2)
    trKkm2 = tracek(i,j,km2)

    !
    ! calc density in sim frame from rest-frame density (rho here)
    !    can show this from: rho_n = T_mu,nu n^mu n^nu = rho_r u_mu u_nu n^mu n^nu
    rho_n = rho * lorentz**2

    ! loop over COMPONENTS of M_i
    do l=1,3

       ! D_j K^j_i term
       term1a = 0._dp
       ! fourth order approx for derivative
       term1a = deriv1fourth(kudip1(1,l),kudip2(1,l),kudim1(1,l),kudim2(1,l),dx) + &
            & deriv1fourth(kudjp1(2,l),kudjp2(2,l),kudjm1(2,l),kudjm2(2,l),dx) + &
            & deriv1fourth(kudkp1(3,l),kudkp2(3,l),kudkm1(3,l),kudkm2(3,l),dx)

       term1b = 0._dp; term1c = 0._dp;
       do m=1,3
          ! loop for sums in christoffel terms
          gamma1b_i = Chrsijk(m,1,m,i,j,k)
          gamma1b_ii = Chrsijk(m,2,m,i,j,k)
          gamma1b_iii = Chrsijk(m,3,m,i,j,k)
          term1b = term1b + ( gamma1b_i * kud(1,l) + gamma1b_ii * kud(2,l) + gamma1b_iii * kud(3,l) )

          gamma1c_i = Chrsijk(1,l,m,i,j,k)
          gamma1c_ii = Chrsijk(2,l,m,i,j,k)
          gamma1c_iii = Chrsijk(3,l,m,i,j,k)
          term1c = term1c + ( gamma1c_i * kud(m,1) + gamma1c_ii * kud(m,2) + gamma1c_iii * kud(m,3) )

          !
          ! calculate K^ij for A^ij for cmo testing
          kup(m,l) = 0._dp
          do n=1,3
            kup(m,l) = kup(m,l) + gup(m,n) * kud(l,n)
          enddo

       enddo
       term1(l) = term1a + term1b - term1c

       ! D_i trK term
       term2(l) = 0._dp
       select case(l)
       case(1)
          ! x-deriv of trK
          term2(l) = deriv1fourth(trKip1,trKip2,trKim1,trKim2,dx)
       case(2)
          ! y-deriv of trK
          term2(l) = deriv1fourth(trKjp1,trKjp2,trKjm1,trKjm2,dx)
       case(3)
          ! z-deriv of trK
          term2(l) = deriv1fourth(trKkp1,trKkp2,trKkm1,trKkm2,dx)
       end select

       ! S_i term
       term3(l) = gdown(l,1) * vel1 + gdown(l,2) * vel2 + gdown(l,3) * vel3
       ! note in the below it's really rho*h, but this is rho * (1 + eps + P/rho) = rho * (1 + K*rho + K) = rho + K*rho^2 + K*rho
       !      so long as K*rho<<rho then this is fine
       term3(l) = 8._dp * pi * term3(l) * rho * lorentz**2 ! this should be rest-frame density (if not, this term is ZERO)
       Mom(l)   = term1(l) - term2(l) - term3(l)

    enddo

    ! calc mag = sqrt{M^i M_i}
    ! take sqares of indivudal terms in Mom to get [M]
    magMom2 = 0._dp; term1mag = 0._dp; term2mag = 0._dp; term3mag = 0._dp

    do l=1,3
       do m=1,3
          magMom2 = magMom2 + gup(l,m) * Mom(l) * Mom(m)
          term1mag = term1mag + gup(l,m) * term1(l) * term1(m)
          term2mag = term2mag + gup(l,m) * term2(l) * term2(m)
          term3mag = term3mag + gup(l,m) * term3(l) * term3(m)
       enddo
    enddo
    Momescale2 =  term1mag + term2mag + term3mag  ! save this for L1 errors

    !
    ! get Hamiltonian constraint violation
    trKijk = tracek(i,j,k)
    call calc_ham_violation(tracer,trKijk,kud,rho_n,kudud,Hamraw,Hamescale2)

    !
    ! Get A_ij and A^ij for cmo testing
    Add = kdown - gdown*trKijk / 3._dp
    AUU = kup - gup*trKijk / 3._dp

    ! calculate rho_Aij -- amplitude of tensor inhomogeneities
    rhoAij = (kudud - (1._dp/3._dp) * trKijk**2) / (16._dp * pi)


  end subroutine constraint_violation_inloop




  !
  ! a routine to calculate the L1 error (raw and/or rel) for the constraint violation
  !       and write these to a file -- called from AFTER the spatial loop in ricci
  !
  subroutine get_L1_constraint_violation(nx,it,time,nargs,Hamraw,Hescale2,magMom2,Momescale2)
    integer, intent(in) :: nx,it,nargs
    real(c_double), intent(in) :: time
    real(c_double), dimension(nx,nx,nx), intent(in) :: Hamraw,Hescale2,magMom2,Momescale2

    real(c_double) :: HamL1raw,HamL1rel,MomL1raw,MomL1rel,tmprad,N
    character(len=100)  :: message,loc
    loc = " get_L1_constraint_violation"

    ! ------------------------------------------------------
    !
    ! calculate L_1 error for Ham and Mom
    !
    ! ------------------------------------------------------

    N = nx**3

    HamL1raw  = sum(abs(Hamraw)) / N !H_L1num
    HamL1rel  = sum(abs(Hamraw) / sqrt(Hescale2)) / N  !H_L1num / H_L1denom

    MomL1raw    = sum(sqrt(magMom2)) / N !Mom_L1num
    !if (Mom_L1denom==0.) then
       ! sometimes, e.g. for FLRW, this is zero and we get NaN (not for Ham)
    !   MomL1rel = Mom_L1num
    !else
    MomL1rel = sum(sqrt(magMom2 / Momescale2)) / N  !Mom_L1num / Mom_L1denom
    !endif

    ! ------------------------------------------------------
    !
    ! write scalar (L1 error) data for H and M
    !
    ! ------------------------------------------------------
    ! for writing; we only need this when using ifort, but it won't affect gfortran
    tmprad = 0._dp
    if (rawHM) then
       !
       ! only print/output RAW violation

       write(message,"(a,ES15.8)") "   --> L1 H violation = ",HamL1raw
       call print_info(message,loc)
       write(message,"(a,ES15.8)") "   --> L1 M violation = ",MomL1raw
       call print_info(message,loc)

       ! last string is "domain_type" -- here this is always all box
       call write_avg((/ HamL1raw /),'Ham_L1error_raw',it,time,nargs,1,tmprad,"all")
       call write_avg((/ MomL1raw /),'Mom_L1error_raw',it,time,nargs,1,tmprad,"all")

    elseif (rawandrelHM) then
       !
       ! print/output RAW + REL violation
       write(message,"(a,ES15.8)") "   --> L1 H violation = ",HamL1raw
       call print_info(message,loc)
       write(message,"(a,ES15.8)") "   --> L1 M violation = ",MomL1raw
       call print_info(message,loc)

       write(message,"(a,ES15.8)") "   --> L1 H/[H] violation = ",HamL1rel
       call print_info(message,loc)
       write(message,"(a,ES15.8)") "   --> L1 M/[M] violation = ",MomL1rel
       call print_info(message,loc)

       call write_avg((/ HamL1rel /),'Ham_L1error_rel',it,time,nargs,1,tmprad,"all")
       call write_avg((/ MomL1rel /),'Mom_L1error_rel',it,time,nargs,1,tmprad,"all")
       call write_avg((/ HamL1raw /),'Ham_L1error_raw',it,time,nargs,1,tmprad,"all")
       call write_avg((/ MomL1raw /),'Mom_L1error_raw',it,time,nargs,1,tmprad,"all")

    else
       !
       ! only print/output REL violation
       write(message,"(a,ES15.8)") "   --> L1 H/[H] violation = ",HamL1rel
       call print_info(message,loc)
       write(message,"(a,ES15.8)") "   --> L1 M/[M] violation = ",MomL1rel
       call print_info(message,loc)

       call write_avg((/ HamL1rel /),'Ham_L1error_rel',it,time,nargs,1,tmprad,"all")
       call write_avg((/ MomL1rel /),'Mom_L1error_rel',it,time,nargs,1,tmprad,"all")

    endif

  end subroutine get_L1_constraint_violation





  !
  ! a routine to calculate the L2 norm (rel ONLY) for the constraint violation
  !       and write these to a file -- called from AFTER the spatial loop in ricci
  !
  subroutine get_L2_constraint_violation(nx,it,time,nargs,Hamraw,Hescale2,magMom2,Momescale2)
    integer, intent(in) :: nx,it,nargs
    real(c_double), intent(in) :: time
    real(c_double), dimension(nx,nx,nx), intent(in) :: Hamraw,Hescale2,magMom2,Momescale2

    real(c_double) :: HamL2rel,MomL2rel,HamL2raw,MomL2raw
    real(c_double) :: tmprad
    character(len=100)  :: message,loc
    loc = " get_L2_constraint_violation"

    ! ------------------------------------------------------
    !
    ! calculate L_2 norm for Ham and Mom
    !
    ! ------------------------------------------------------


    HamL2raw   = sqrt(sum( Hamraw**2 ))
    HamL2rel  =  sqrt( sum( Hamraw**2/Hescale2 ) )    !HamL2raw / H_L2denom

    MomL2raw   = sqrt(sum(magMom2))
    if (sum(Momescale2)==0.) then
       ! sometimes, e.g. for FLRW, this is zero and we get NaN (not for Ham)
       MomL2rel = MomL2raw
    else
       MomL2rel = sqrt( sum( magMom2/Momescale2 ) )    !MomL2raw / Mom_L2denom
    endif

    ! ------------------------------------------------------
    !
    ! write scalar (L2 norm) data for H and M
    !
    ! ------------------------------------------------------
    ! a workaround for intel compilers, don't like passing 0. to real(c_double)
    tmprad = 0._dp
    if (rawHM) then
       !
       ! only print/output RAW violation

       write(message,"(a,ES15.8)") "   --> L2 H violation = ",HamL2raw
       call print_info(message,loc)
       write(message,"(a,ES15.8)") "   --> L2 M violation = ",MomL2raw
       call print_info(message,loc)

       call write_avg((/ HamL2raw /),'Ham_L2error_raw',it,time,nargs,1,tmprad,"all")
       call write_avg((/ MomL2raw /),'Mom_L2error_raw',it,time,nargs,1,tmprad,"all")

    elseif (rawandrelHM) then
       !
       ! print/output RAW + REL violation
       write(message,"(a,ES15.8)") "   --> L2 H violation = ",HamL2raw
       call print_info(message,loc)
       write(message,"(a,ES15.8)") "   --> L2 M violation = ",MomL2raw
       call print_info(message,loc)

       write(message,"(a,ES15.8)") "   --> L2 H/[H] violation = ",HamL2rel
       call print_info(message,loc)
       write(message,"(a,ES15.8)") "   --> L2 M/[M] violation = ",MomL2rel
       call print_info(message,loc)

       call write_avg((/ HamL2rel /),'Ham_L2error_rel',it,time,nargs,1,tmprad,"all")
       call write_avg((/ MomL2rel /),'Mom_L2error_rel',it,time,nargs,1,tmprad,"all")
       call write_avg((/ HamL2raw /),'Ham_L2error_raw',it,time,nargs,1,tmprad,"all")
       call write_avg((/ MomL2raw /),'Mom_L2error_raw',it,time,nargs,1,tmprad,"all")

    else
       !
       ! only print/output REL violation
       write(message,"(a,ES15.8)") "   --> L2 H/[H] violation = ",HamL2rel
       call print_info(message,loc)
       write(message,"(a,ES15.8)") "   --> L2 M/[M] violation = ",MomL2rel
       call print_info(message,loc)

       call write_avg((/ HamL2rel /),'Ham_L2error_rel',it,time,nargs,1,tmprad,"all")
       call write_avg((/ MomL2rel /),'Mom_L2error_rel',it,time,nargs,1,tmprad,"all")

    endif

  end subroutine get_L2_constraint_violation



  ! ======================================================================================================
  !
  ! Calculate the constraint violation in the FLUID FRAME
  !
  ! see Buchert 2000, eqs 4 for constraints in fluid decomposition
  !      - these are in terms of covD of quantities where covD is w.r.t the g_ij of fluid HS's
  !      - we can project this into our sim frame: D_i --> h_i^j D_j, and corrections are O(v^2) for g_ij = O(1)
  !      - so for now this is just the D_i in sim frame, can add correction terms later
  !
  ! Ham constraint: general R passed in, can be 3R or fluidR depending on tdeirvs=False/True
  !                 assumed zero vorticity
  !
  ! ======================================================================================================
  subroutine get_fluid_frame_constraints(nx,dx,gij,sigmaUd,theta,sigma2,ricciscalar,rho,Gammas,FluidHam,&
      & FluidHamescale2,magFluidMom2,FluidMomescale2)
      integer, intent(in) :: nx
      real(c_double), intent(in) :: dx,gij(6,nx,nx,nx),sigmaUd(4,4,nx,nx,nx),theta(nx,nx,nx)
      real(c_double), intent(in) :: sigma2(nx,nx,nx),ricciscalar(nx,nx,nx),rho(nx,nx,nx)
      real(c_double), intent(in) :: Gammas(3,3,3,nx,nx,nx) ! Christoffels

      real(c_double), intent(out) :: magFluidMom2(nx,nx,nx),FluidMomescale2(nx,nx,nx),&
        & FluidHam(nx,nx,nx),FluidHamescale2(nx,nx,nx)

      real(c_double) :: gdown(3,3),FMom(3),magMom,term1(3),term2(3),term1mag,term2mag
      real(c_double) :: disigUd(3,3,3),ditheta(3),dxtheta,dytheta,dztheta,gamijk(3,3,3)
      real(c_double) :: sigUdijk(3,3),sigUdip1(3,3),sigUdip2(3,3),sigUdim1(3,3),sigUdim2(3,3)
      real(c_double) :: sigUdjp1(3,3),sigUdjp2(3,3),sigUdjm1(3,3),sigUdjm2(3,3)
      real(c_double) :: sigUdkp1(3,3),sigUdkp2(3,3),sigUdkm1(3,3),sigUdkm2(3,3)

      integer :: i,j,k,m,ii,jj
      integer :: ip1,jp1,kp1,im1,jm1,km1,ip2,im2,jp2,jm2,kp2,km2 ! for periodic boundaries
      character(len=100) :: loc
      loc = "  get_fluid_frame_constraints"

      !$omp parallel do default (none) &
      !$omp shared(gij,magFluidMom2,FluidMomescale2,FluidHam,FluidHamescale2) &
      !$omp shared(sigmaUd,theta,Gammas,rho,sigma2,ricciscalar,nx,dx) &
      !$omp private(ip1,jp1,kp1,im1,jm1,km1,ip2,im2,jp2,jm2,kp2,km2,gdown,FMom) &
      !$omp private(m,ii,jj) &
      !$omp private(ditheta,dxtheta,dytheta,dztheta,gamijk,magMom) &
      !$omp private(sigUdijk,sigUdip1,sigUdip2,sigUdim1,sigUdim2,term1,term2,term1mag,term2mag) &
      !$omp private(sigUdjp1,sigUdjp2,sigUdjm1,sigUdjm2,sigUdkp1,sigUdkp2,sigUdkm1,sigUdkm2) &
      !$omp private(disigUd)
      do k=1,nx
          do j=1,nx
              do i=1,nx

                  ! --------------------------------------------------------
                  !
                  ! Get stencils for spatial deriv of sigmaUd and theta
                  !
                  ! --------------------------------------------------------
                  call apply_periodic(i,ip1,im1,nx)
                  call apply_periodic(j,jp1,jm1,nx)
                  call apply_periodic(k,kp1,km1,nx)
                  call apply_periodic_fourth(i,ip2,im2,nx)
                  call apply_periodic_fourth(j,jp2,jm2,nx)
                  call apply_periodic_fourth(k,kp2,km2,nx)

                  !
                  ! calculate spatial derivatives of shear and theta terms
                  !
                  dxtheta = deriv1fourth(theta(ip1,j,k),theta(ip2,j,k),theta(im1,j,k),theta(im2,j,k),dx)
                  dytheta = deriv1fourth(theta(i,jp1,k),theta(i,jp2,k),theta(i,jm1,k),theta(i,jm2,k),dx)
                  dztheta = deriv1fourth(theta(i,j,kp1),theta(i,j,kp2),theta(i,j,km1),theta(i,j,km2),dx)
                  ditheta = (/ dxtheta, dytheta, dztheta /)

                  !
                  ! We need all derivatives of all 3x3 components of sigUd
                  sigUdijk = sigmaUd(2:4,2:4,i,j,k)

                  sigUdip1 = sigmaUd(2:4,2:4,ip1,j,k)
                  sigUdip2 = sigmaUd(2:4,2:4,ip2,j,k)
                  sigUdim1 = sigmaUd(2:4,2:4,im1,j,k)
                  sigUdim2 = sigmaUd(2:4,2:4,im2,j,k)

                  sigUdjp1 = sigmaUd(2:4,2:4,i,jp1,k)
                  sigUdjp2 = sigmaUd(2:4,2:4,i,jp2,k)
                  sigUdjm1 = sigmaUd(2:4,2:4,i,jm1,k)
                  sigUdjm2 = sigmaUd(2:4,2:4,i,jm2,k)

                  sigUdkp1 = sigmaUd(2:4,2:4,i,j,kp1)
                  sigUdkp2 = sigmaUd(2:4,2:4,i,j,kp2)
                  sigUdkm1 = sigmaUd(2:4,2:4,i,j,km1)
                  sigUdkm2 = sigmaUd(2:4,2:4,i,j,km2)

                  ! loop over components and calculate x,y,z derivs and store
                  do ii=1,3
                      do jj=1,3
                          ! x-derivative
                          disigUd(1,ii,jj) = deriv1fourth(sigUdip1(ii,jj),sigUdip2(ii,jj),sigUdim1(ii,jj),sigUdim2(ii,jj),dx)

                          ! y-derivative
                          disigUd(2,ii,jj) = deriv1fourth(sigUdjp1(ii,jj),sigUdjp2(ii,jj),sigUdjm1(ii,jj),sigUdjm2(ii,jj),dx)

                          ! z-derivative
                          disigUd(3,ii,jj) = deriv1fourth(sigUdkp1(ii,jj),sigUdkp2(ii,jj),sigUdkm1(ii,jj),sigUdkm2(ii,jj),dx)
                      enddo
                  enddo

                  !
                  ! build the violation
                  !
                  gamijk = Gammas(:,:,:,i,j,k)
                  term1 = 0._dp; term2 = 0._dp
                  FMom  = 0._dp
                  ! violation component loop
                  do m=1,3
                      ! -2/3 \pd_m \theta
                      FMom(m) = - (2._dp / 3._dp) * ditheta(m)
                      ! store the term for the magnitude calculation
                      term2(m) = (2._dp / 3._dp) * ditheta(m)

                      ! first sum loop
                      do ii=1,3
                          ! \pd_i sig^i_m
                          FMom(m) = FMom(m) + disigUd(ii,ii,m)

                          ! store the term for the magnitude calculation
                          term1(m) = term1(m) + disigUd(ii,ii,m)

                          ! second sum loop
                          do jj=1,3
                              ! \Gam^i_{ik} \sig^k_m - \Gam^k_{im} \sig^i_k
                              FMom(m) = FMom(m) + gamijk(ii,ii,jj) * sigUdijk(jj,m) - gamijk(jj,ii,m) * sigUdijk(ii,jj)

                              ! store the term for the magnitude calculation
                              term1(m) = term1(m) + gamijk(ii,ii,jj) * sigUdijk(jj,m) - gamijk(jj,ii,m) * sigUdijk(ii,jj)
                          enddo

                      enddo

                  enddo
                  !
                  ! get the magnitude of the violation
                  !
                  call get_metric_at_pos(i,j,k,nx,gij,gdown)
                  magMom = 0._dp; term1mag = 0._dp; term2mag = 0._dp
                  do ii=1,3
                      do jj=1,3
                          magMom = magMom + gdown(ii,jj) * FMom(ii) * FMom(jj)

                          term1mag = term1mag + gdown(ii,jj) * term1(ii) * term1(jj)
                          term2mag = term2mag + gdown(ii,jj) * term2(ii) * term2(jj)
                      enddo
                  enddo
                  magFluidMom2(i,j,k)    = magMom
                  FluidMomescale2(i,j,k) = term1mag + term2mag

                  !
                  ! Get the hamiltonian constraint here and energy scale
                  !
                  FluidHam(i,j,k)        = (2._dp / 3._dp) * theta(i,j,k)**2 - 2._dp * sigma2(i,j,k) &
                    & + ricciscalar(i,j,k) - 16._dp*pi*rho(i,j,k)
                  FluidHamescale2(i,j,k) = ( (2._dp / 3._dp) * theta(i,j,k)**2 )**2 + (2._dp * sigma2(i,j,k) )**2 &
                    & + ricciscalar(i,j,k)**2 + ( 16._dp*pi*rho(i,j,k) )**2

              enddo
          enddo
      enddo
      !$omp end parallel do

  end subroutine get_fluid_frame_constraints




end module violation
