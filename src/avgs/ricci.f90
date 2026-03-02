!-------------------------------------------------------------------------
!
!  Module to compute the ricci tensor (and other things) from the metric
!
!     This module implements the new averaging scheme for general spacetime
!     foliations from Buchert, Mourier, & Roy (2019)
!
!-------------------------------------------------------------------------
module ricci
    use options, only:eps_t,dp,c_double,clen,tinitial,dit,dtfac,nord_dt,tderivs,&
        & constraints_l1norm,writeomegas,writeaD,writehubble,writedelta,writetilde,&
        & write3D,shear_3Dout,Rdd_3Dout,Aij_3Dout,EUd_3Dout,shearUd_3Dout,gdd_3Dout,&
        & fluid_constraints_3Dout,domain_type
    use roots, only:get_ricci_component,get_backreaction_omegas,&
        & get_expansion_shear_vort,get_fluid_curvature,get_christoffels_trK_allgrid_alltime
    use weyl, only:get_electric_weyl
    use spines, only:calc_average_within_radius,calc_average_within_radius_manyspheres
    use gauge, only:get_dtalp
    use manipulations, only:get_metric_at_pos,trace,inv3x3,get_metric_at_stencil,&
        & calc_up_down_ND
    use tripreports, only:write_avg,write_hdf5
    use violation, only:constraint_violation_inloop,get_L1_constraint_violation,&
        & get_fluid_frame_constraints
    use init, only:initialise
    use prints, only:print_info,print_error
    implicit none

contains
  !
  ! This is the main calculation routine that calls other routines to calculate:
  !      - Ricci tensor and scalar (4-D and 3-D)
  !      - Christoffel symbols
  !      - trK, expansion, shear, vorticity, acceleration
  !      - perform spatial averaging, cosmological parameters and backreaction
  !      - constraint violation
  !      - write main 3D output files
  !
  subroutine compute_ricci(nx,nt,dx,xmax,xmin,rad,nspheres,gij,kij,rho,vel0,vel1,vel2,alp,times,&
       gotrho,gotalp,dens)
    integer, intent(in) :: nx,nt,nspheres
    ! NOTE: ordering of times is e.g. (t_n-3, t_n-2, t_n-1, t_n) for nord_dt = nt = 4
    real(c_double), intent(in) :: gij(6,nx,nx,nx,nt),kij(6,nx,nx,nx,nt),alp(nx,nx,nx)
    real(c_double), intent(in) :: times(nt),dx,rad,xmax,xmin,dens(nx,nx,nx)
    real(c_double), intent(in), dimension(nx,nx,nx,nt) :: vel0,vel1,vel2
    real(c_double), intent(inout) :: rho(nx,nx,nx)     ! this is (inout) because we may need to set it if gotrho=False
    logical, intent(in) :: gotalp,gotrho               ! whether we have hdf5 primrho (lapse) data, if false then calculate it from dens (dont do thing that require alp)

    real(c_double) :: xvals(nx),dt,xi,xj,xk,time
    real(c_double) :: lorentzijk,tracekijk,tracerfluidijk,tracerijk,Rvsqijk
    real(c_double) :: thetaijk,sigma2ijk,detgijk,rhoijk,w2ijk,alpijk
    real(c_double) :: avgrho_allbox,avgrho_allbox_notilde,avglorentz_allbox,hub_allbox_notilde
    real(c_double) :: fac,fac2,thetaavg,traceravg,rhoavg,theta2avg,sigma2avg,vort2avg,tmprad
    real(c_double), dimension(3)           :: dialp,dtudi,diudt
    real(c_double), dimension(4)           :: adijk,umudijk,aUijk,umuUijk
    real(c_double), dimension(3,3)         :: gdown,gup,threeRij,Addijk,AUUijk
    real(c_double), dimension(4,4)         :: sigUUijk,wddijk,sigUdijk,sigddijk
    real(c_double), dimension(4,4)         :: fourRddijk,Elecddijk,ElecUdijk,gUUmunu
    real(c_double), dimension(nspheres)    :: hub,Qd,omegaQ
    real(c_double), dimension(nspheres)    :: omegaR,omegam,aDb,delta
    real(c_double), dimension(3,nspheres)  :: randorigins
    real(c_double), dimension(13,3,3)      :: gdown_atstencil,kdown_atstencil,kud_atstencil
    real(c_double) :: fvelU_atstencil(13,3),dtudt,dtalp
    real(c_double) :: rxx,rxy,rxz,ryy,ryz,rzz,rtensor(6),fluidRijk,fourRicci
    real(c_double) :: magaijk,rhoAij_ijk,Esqijk

    integer, parameter :: navgs = 7
    real(c_double), dimension(navgs,nspheres) :: avgs
    real(c_double), dimension(:,:,:), allocatable :: detg,theta,sigma2,w2,lorentz,maga,rhoAij
    real(c_double), dimension(:,:,:), allocatable :: tracerfluid,tracer,tracefourr,Esq
    real(c_double), dimension(:,:,:), allocatable :: magMom2,Momescale2,Hamraw,Hamescale2
    real(c_double), dimension(:,:,:), allocatable :: magFluidMom2,FluidMomescale2,FluidHam,FluidHamescale2
    real(c_double), dimension(:,:,:,:),     allocatable :: datas,tracek,kij_tn,gij_tn,acceld,fveld
    real(c_double), dimension(:,:,:,:,:),   allocatable :: sigmadd,fourRdd,Aijdd,ElecUd,sigmaUd
    real(c_double), dimension(:,:,:,:,:,:),   allocatable :: gamijk

    character(len=clen), allocatable :: descriptors(:)
    character(len=100) :: message,loc

    integer :: k,j,i,it,ndat,ni
    logical :: notearly  ! if True, then we calc \Theta. if False, it is too early (i.e. we can't take time derivs yet)
    loc      = "Ricci"   ! routine location for info/errors
    notearly = .False.

    ! -----------------------------------------------------
    !
    ! Allocate memory for default HDF5 output and for things we need (e.g. trK, Christoffels)
    !
    ! -----------------------------------------------------
    allocate(tracerfluid(nx,nx,nx),lorentz(nx,nx,nx),sigma2(nx,nx,nx),&
         & theta(nx,nx,nx),detg(nx,nx,nx),magMom2(nx,nx,nx),Momescale2(nx,nx,nx),&
         & Hamraw(nx,nx,nx),Hamescale2(nx,nx,nx),tracer(nx,nx,nx),rhoAij(nx,nx,nx),&
         & gamijk(3,3,3,nx,nx,nx),tracek(nx,nx,nx,nt),w2(nx,nx,nx),&
         & tracefourr(nx,nx,nx),maga(nx,nx,nx))
    allocate(kij_tn(6,nx,nx,nx),gij_tn(6,nx,nx,nx)) ! K_{ij}, g_{ij} at the current time

    ! -----------------------------------------------------
    !
    ! Allocate some memory for additional output
    !
    ! -----------------------------------------------------

    ! Extra memory for all components of shear tensor output
    if (shear_3Dout) then
        if (shearUd_3Dout) then
            ! We want sigma^{mu}_{nu}
            allocate(sigmaUd(4,4,nx,nx,nx))
        else
            ! We want sigma_{mu,nu}
            allocate(sigmadd(4,4,nx,nx,nx))
        endif
    endif
    ! Extra memory for 4-Ricci tensor components
    if (Rdd_3Dout)   allocate(fourRdd(4,4,nx,nx,nx))
    ! Extra memory for A_ij 3D output
    if (Aij_3Dout)   allocate(Aijdd(3,3,nx,nx,nx))
    ! Extra memory for E^{mu}_{nu} 3D output and E^2 scalar
    if (EUd_3Dout)   allocate(ElecUd(4,4,nx,nx,nx),Esq(nx,nx,nx))

    ! -----------------------------------------------------
    !
    ! Intialise some grid, sphere, things etc. and filenames
    !
    ! -----------------------------------------------------
    time = times(nt) ! the current time
    call initialise(nx,nspheres,rad,time,xmin,xmax,dx,gotrho,xvals,randorigins,dt,&
         & it,notearly)

    !
    ! initialise our averages
    avgrho_allbox = 0._dp; avgrho_allbox_notilde = 0._dp; avglorentz_allbox = 0._dp
    hub_allbox_notilde = 0._dp
    avgs = 0._dp; omegam = 0._dp; omegaR = 0._dp; omegaQ = 0._dp

    !
    ! A workaround for intel compilers for passing hard-coded zero radius for avgs that are always all-box
    tmprad = 0._dp

    ! --------------------------------------------------------------------------------------------
    !                                    user information
    ! --------------------------------------------------------------------------------------------
    !
    if (tderivs) then
       if (notearly) then
          write(message,"(a,i1)") "Taking time derivatives of order: ",nord_dt
          call print_info(message,loc)
       else
          call print_info("Too early to take backward time derivatives. Using linear approx. for expansion & shear (for now) ",loc)
       endif
    else
       call print_info("Not taking time derivatives. Using linear approx. for expansion & shear instead ",loc)
       ! Check we're not trying to calculate Weyl with tderivs=False; code will fail if so
       if (EUd_3Dout) then
           call print_error("You need to set tderivs=True to calculate the Weyl tensor. Re-set and re-run.",2,loc)
       endif
    endif
    call print_info("Using avergaing formalism of Buchert+(2019) ",loc)
    !
    ! --------------------------------------------------------------------------------------------

    !
    ! Get some things we need to get before the main loop
    !
    call print_info("  Calculating Christoffel symbols and trace(K_ij) ... ",loc)
    kij_tn = kij(:,:,:,:,nt) ! to save accessing memory in this loop
    gij_tn = gij(:,:,:,:,nt)
    ! 2nd instance of 'nt' in below call is 'ti', i.e. index in time dim of 'current time'
    call get_christoffels_trK_allgrid_alltime(nx,nt,nt,dx,gij,kij,gamijk,tracek)

    call print_info("  Calculating 3R, expansion, shear, constraint violation, and averaging ... ",loc)
    !$omp parallel do default (none) &
    !$omp shared(nx,nt,dx,time,times,it,tracerfluid,rho,dens,tracek,kij,kij_tn,gij_tn) &
    !$omp shared(notearly,gotrho,gotalp,vel0,vel1,vel2,lorentz,xvals) &
    !$omp shared(sigma2,w2,dt,gamijk,magMom2,Momescale2,Hamraw,Hamescale2,tmprad) &
    !$omp shared(nspheres,rad,randorigins,detg,alp,theta,tracer,rhoAij,sigmadd,fourRdd) &
    !$omp shared(Aijdd,maga,xmax,tracefourr,acceld,fveld,sigmaUd,ElecUd,Esq) &
    !$omp private(i,j,k,gdown,gup,dtalp,gdown_atstencil,kdown_atstencil,Addijk,AUUijk) &
    !$omp private(rxx,rxy,rxz,ryy,ryz,rzz,kud_atstencil,fourRicci,fourRddijk,gUUmunu) &
    !$omp private(aUijk,adijk,magaijk,Rvsqijk,umudijk,Elecddijk,ElecUdijk,sigUdijk) &
    !$omp private(fvelU_atstencil,dialp,dtudt,dtudi,diudt,fluidRijk,rtensor,threeRij) &
    !$omp private(detgijk,thetaijk,sigma2ijk,sigddijk,sigUUijk,w2ijk,lorentzijk,rhoijk,rhoAij_ijk) &
    !$omp private(alpijk,fac,fac2,thetaavg,traceravg,rhoavg,tracekijk,tracerfluidijk,tracerijk) &
    !$omp private(theta2avg,sigma2avg,xi,xj,xk,vort2avg,wddijk,umuUijk,Esqijk) &
    !$omp reduction(+:avgs,avgrho_allbox,avgrho_allbox_notilde,hub_allbox_notilde,avglorentz_allbox)
    do k=1,nx
       do j=1,nx
          do i=1,nx
             !
             ! reduce need to access memory throughout loop
             xi = xvals(i); xj = xvals(j); xk = xvals(k)
             !
             ! Get some metric and K_ij at *spatial* finite differencing stencil, and g^{ij}
             !
             call get_metric_at_stencil(i,j,k,nx,gij_tn,gdown_atstencil)
             gdown = gdown_atstencil(1,:,:)
             call inv3x3(gdown,gup,detgijk)
             call get_metric_at_stencil(i,j,k,nx,kij_tn,kdown_atstencil)

             !
             ! Get the relevant scalars to average
             !
             ! 1. get lapse time derivatives
             ! 2. get fluid expansion \Theta
             ! 3. get shear \sigma^2
             !      -- note 2. and 3. are correct to O(v^2) if notearly=False or tderivs=False in options
             !
             alpijk    = alp(i,j,k)
             tracekijk = tracek(i,j,k,nt)
             call get_dtalp(time,alpijk,tracekijk,dtalp)

             call get_expansion_shear_vort(i,j,k,nx,nt,dx,dit*dt,notearly,gdown_atstencil,kdown_atstencil(1,:,:),&
                & gamijk(:,:,:,i,j,k),alp,dtalp,gup,vel0,vel1,vel2,tracekijk,fvelU_atstencil,umudijk,dialp,&
                & dtudt,dtudi,diudt,sigma2ijk,thetaijk,w2ijk,lorentzijk,sigddijk,sigUUijk,wddijk)
             if (shear_3Dout) then
                 if (shearUd_3Dout) then
                     ! Raise an index of sigma_{mu,nu} -- lots of loops but unavoidable while staying general
                     !    first build g^{mu,nu} to pass in
                     gUUmunu      = 0._dp
                     gUUmunu(1,1) = - 1._dp / alpijk**2      ! g^{00} = -1/alp^2
                     gUUmunu(2:4,2:4) = gup                  ! g^{ij} = \gamma^{ij} (for beta^i=0)
                     call calc_up_down_ND(4,sigddijk,gUUmunu,sigUdijk)
                 endif
             endif

             if (gotrho) then
                !
                ! we have primrho data - store this position's data
                rhoijk = rho(i,j,k)
             else
                !
                ! we have no primrho data: calculate it using conserved dens (D - see GRHydro doc)
                rhoijk = dens(i,j,k) / ( sqrt(detgijk) * lorentzijk )
             endif
             !
             !     Get curvature scalar in fluid rest frame \mathcal{R}
             !
             !
             ! First, get 3-Ricci curvature + put into (6) for trace and (3,3) for fluidR call
             !
             call get_ricci_component(i,j,k,nx,dx,gij_tn,1,1,rxx,gamijk)
             call get_ricci_component(i,j,k,nx,dx,gij_tn,1,2,rxy,gamijk)
             call get_ricci_component(i,j,k,nx,dx,gij_tn,1,3,rxz,gamijk)
             call get_ricci_component(i,j,k,nx,dx,gij_tn,2,2,ryy,gamijk)
             call get_ricci_component(i,j,k,nx,dx,gij_tn,2,3,ryz,gamijk)
             call get_ricci_component(i,j,k,nx,dx,gij_tn,3,3,rzz,gamijk)
             rtensor   = (/ rxx, rxy, rxz, ryy, ryz, rzz /)
             tracerijk = trace(rtensor,gdown)
             threeRij(1,1) = rxx; threeRij(2,2) = ryy; threeRij(3,3) = rzz
             threeRij(1,2) = rxy; threeRij(2,1) = rxy
             threeRij(1,3) = rxz; threeRij(3,1) = rxz
             threeRij(2,3) = ryz; threeRij(3,2) = ryz
             !
             if (tderivs .and. notearly) then
                call get_fluid_curvature(nx,nt,dx,dit*dt,i,j,k,gdown_atstencil,gup,kij,kdown_atstencil,&
                     & tracek,alp,threeRij,gamijk(:,:,:,i,j,k),thetaijk,lorentzijk,fvelU_atstencil,dtudt,dtudi,&
                     & dtalp,dialp,diudt,kud_atstencil,aUijk,adijk,magaijk,fluidRijk,fourRicci,fourRddijk,Rvsqijk)
                !fluidRijk = fluid_restframe_curvature(rhoijk,thetaijk,sigma2ijk,w2ijk)
             else
                fourRicci = 0._dp; fluidRijk = 0._dp; fourRddijk = 0._dp
                aUijk = 0._dp; adijk = 0._dp; magaijk = 0._dp
             endif
             !
             ! Get Electric Weyl if we want to
             if (EUd_3Dout) then
                 ! build u^mu = (u^0, u^i) = (Gam/alp, u^i)
                 umuUijk = (/ lorentzijk/alpijk, fvelU_atstencil(1,1), fvelU_atstencil(1,2), fvelU_atstencil(1,3) /)
                 ! second instance of 'nt' in below call is 'ti'
                 call get_electric_weyl(nx,nt,nt,i,j,k,dx,dit*dt,gij_tn,kij,alp,gamijk,gdown,&
                    & gup,kdown_atstencil(1,:,:),fourRicci,fourRddijk,umuUijk,umudijk,Elecddijk,Esqijk)
                 ! Now raise one index -- this is too many loops but I don't think we can avoid it
                 !    first build g^{mu,nu} to pass in
                 gUUmunu      = 0._dp
                 gUUmunu(1,1) = - 1._dp / alpijk**2      ! g^{00} = -1/alp^2
                 gUUmunu(2:4,2:4) = gup                  ! g^{ij} = \gamma^{ij} (for beta^i=0)
                 call calc_up_down_ND(4,Elecddijk,gUUmunu,ElecUdijk)
             endif

             !
             ! Build scalars to average according to new Buchert+ (2019) formalism
             !   (note: these are all multiplied by \gamma to convert to <\psi>^b_D from <\psi>_D
             !
             fac       = alpijk / lorentzijk
             fac2      = fac**2
             thetaavg  = lorentzijk * fac * thetaijk
             if (tderivs .and. notearly) then
                 ! We have calculated fluid curvature above; use that
                 traceravg = lorentzijk * fac2 * fluidRijk
             else
                 ! We don't have fluid curvature; use 3Ricci
                 traceravg = lorentzijk * fac2 * tracerijk
             endif
             rhoavg    = lorentzijk * fac2 * rhoijk
             sigma2avg = lorentzijk * fac2 * sigma2ijk
             theta2avg = lorentzijk * fac2 * thetaijk**2
             vort2avg  = lorentzijk * fac2 * w2ijk
             !
             ! Calculate averages within spheres / over whole domain
             !   avgs: theta, R, rho, sigma2, theta2, w2, lorentz (to convert <>^b_D from <>_D)
             !
             call calc_average_within_radius_manyspheres(xi,xj,xk,rad,nspheres,randorigins,&
                  & detgijk,dx,navgs,avgs,(/ thetaavg, traceravg, rhoavg, sigma2avg, theta2avg, &
                  & vort2avg, lorentzijk /))
             !
             ! calculate average density over the whole box; for density contrast \delta
             ! all box avg ==> rad=0., nspheres=1, randorigins=0.
             !
             ! --> include avg rho WITHOUT scaling, for dL calc
             !
             ! and calculate V^b_D = \sum_D \Gamma for global avg density, delta
             !
             call calc_average_within_radius(xi,xj,xk,(/tmprad,tmprad,tmprad/),tmprad,detgijk,dx,&
                  & avgrho_allbox,rhoavg)
             call calc_average_within_radius(xi,xj,xk,(/tmprad,tmprad,tmprad/),tmprad,detgijk,dx,&
                  & avgrho_allbox_notilde,lorentzijk*rhoijk)
             call calc_average_within_radius(xi,xj,xk,(/tmprad,tmprad,tmprad/),tmprad,detgijk,dx,&
                  & avglorentz_allbox,lorentzijk)
             !
             ! Also calculate H_all = <Theta>/3 without tilde (this is useful for Asta's formalism)
             call calc_average_within_radius(xi,xj,xk,(/tmprad,tmprad,tmprad/),tmprad,detgijk,dx,&
                & hub_allbox_notilde,lorentzijk*thetaijk)

             !
             ! Calculate constraint violation
             !
             call constraint_violation_inloop(i,j,k,nx,dx,gdown_atstencil,kdown_atstencil,gamijk,&
                  & tracek(:,:,:,nt),rhoijk,vel0(i,j,k,nt),vel1(i,j,k,nt),vel2(i,j,k,nt),lorentzijk,&
                  & tracerijk,rhoAij_ijk,Addijk,AUUijk,magMom2(i,j,k),Momescale2(i,j,k),Hamraw(i,j,k),&
                  & Hamescale2(i,j,k))

             !
             ! Store some things for default 3D output
             !
             lorentz(i,j,k)     = lorentzijk ! this used to be vsq output -- but this might save some calcs?
             detg(i,j,k)        = detgijk
             tracefourr(i,j,k)  = fourRicci
             tracer(i,j,k)      = tracerijk
             rhoAij(i,j,k)      = rhoAij_ijk
             maga(i,j,k)        = magaijk

             if (writetilde) then
                ! Write tilde 3D quantities
                if (gotrho .eqv. .False.) rho(i,j,k) = fac2 * rhoijk
                theta(i,j,k)       = fac * thetaijk
                sigma2(i,j,k)      = fac2 * sigma2ijk
                w2(i,j,k)          = fac2 * w2ijk
                tracerfluid(i,j,k) = fac2 * fluidRijk !tracerfluidijk
             else
                ! Write regular quantities
                if (gotrho .eqv. .False.) rho(i,j,k) = rhoijk
                theta(i,j,k)       = thetaijk
                sigma2(i,j,k)      = sigma2ijk
                w2(i,j,k)          = w2ijk
                tracerfluid(i,j,k) = fluidRijk !tracerfluidijk
             endif

             !
             ! Store extra 3D output if we want it
             !
             if (shear_3Dout) then
                 if (shearUd_3Dout) then
                     ! We want sigma^mu_nu
                     sigmaUd(:,:,i,j,k) = sigUdijk
                 else
                     ! We want sigma_{mu,nu}
                     sigmadd(:,:,i,j,k) = sigddijk
                 endif
             endif
             if (Rdd_3Dout)   fourRdd(:,:,i,j,k) = fourRddijk
             if (Aij_3Dout)     Aijdd(:,:,i,j,k) = Addijk
             if (EUd_3Dout) then
                 ElecUd(:,:,i,j,k) = ElecUdijk
                 Esq(i,j,k)        = Esqijk
             endif

          enddo
       enddo
    enddo
    !$omp end parallel do

    ! -------------------------------------------------------------------------------
    !
    ! Calculate the constraint violation in the fluid frame
    !      -- we need theta and sigma^i_j at all grid points for this
    !      -- this is why it's an extra loop here, since we won't need it very often
    !
    ! -------------------------------------------------------------------------------
    if (fluid_constraints_3Dout) then
        if (shear_3Dout) then
            if (shearUd_3Dout) then
                ! Allocate the arrays
                allocate(magFluidMom2(nx,nx,nx),FluidMomescale2(nx,nx,nx),&
                    & FluidHam(nx,nx,nx),FluidHamescale2(nx,nx,nx))
                if (tderivs) then
                    ! we have fluid-frame R, pass this
                    call get_fluid_frame_constraints(nx,dx,gij,sigmaUd,theta,sigma2,tracerfluid,rho,gamijk,FluidHam,&
                        & FluidHamescale2,magFluidMom2,FluidMomescale2)
                else
                    ! we only have 3R, pass this
                    call get_fluid_frame_constraints(nx,dx,gij,sigmaUd,theta,sigma2,tracer,rho,gamijk,FluidHam,&
                        & FluidHamescale2,magFluidMom2,FluidMomescale2)
                endif
                !
                ! Also get the L1 error for these
                !
                !if (constraints_l1norm) call get_L1_constraint_violation(nx,it,time,FluidHam,&
                !    & FluidHamescale2,magFluidMom2,FluidMomescale2)
            else
                ! no shearUd_3Dout:
                call print_error("Set shear_3Dout AND shearUd_3Dout to TRUE for fluid constraints",2,loc)
            endif
        else
            ! no shear_3Dout:
            call print_error("Set shear_3Dout AND shearUd_3Dout to TRUE for fluid constraints",2,loc)
        endif
    endif

    ! --------------------------------------------------------------------------------
    !
    !    Volume of comoving domain is:
    !
    !    V^b_D = V_D * <lorentz>_D
    !          = \int_D lorentz \sqrt{gamma} d^3X (1/V_D cancels -- see "user guide")
    !
    ! --------------------------------------------------------------------------------

    call print_info("  Calculating backreaction terms and cosmological parameters ... ",loc)
    call get_backreaction_omegas(it,time,rad,nspheres,randorigins,avgs,navgs,avgrho_allbox,&
         & avgrho_allbox_notilde,avglorentz_allbox,Qd,hub,omegam,omegaR,omegaQ,delta,aDb)
    if (nspheres<10) then
       write(message,"(a,ES14.7)") "      <rho>_D(all) = ",avgrho_allbox
       call print_info(message,loc)
       write(message,"(a,ES14.7)") "          V_D(all) = ",avglorentz_allbox
       call print_info(message,loc)
    endif

    ! --------------------------------------------------------------------------------
    !
    ! Get L1/L2 error for constraints -- and writes to files (scalars)
    !
    ! --------------------------------------------------------------------------------
    if (constraints_l1norm) call get_L1_constraint_violation(nx,it,time,Hamraw,Hamescale2,magMom2,Momescale2)
    !if (constraints_l2norm) call get_L2_constraint_violation(nx,it,time,Hamraw,Hamescale2,magMom2,Momescale2)

    ! --------------------------------------------------------------------------------
    !
    !  ----- WRITE 3D DATA -----
    !
    ! --------------------------------------------------------------------------------
    if (write3D) then
        call print_info(" Writing 3D (HDF5) data ... ",loc)

        ndat = 19                           ! number of default 3D data arrays we want to write
        if (shear_3Dout) then
            if (shear_3Dout) then
                ! output sigma^{mu}_{nu} -- NOT symmetric so need all x16
                ndat = ndat + 16
            else
                ! outputting sigma_{mu,nu} -- symmetric
                ndat = ndat + 10
            endif
        endif
        if (Rdd_3Dout)     ndat = ndat + 10 ! add 10x independent components of 4R_mu,nu
        if (Aij_3Dout)     ndat = ndat + 6  ! add 6x independent components of A_ij
        if (EUd_3Dout)     ndat = ndat + 17 ! add all 16x compontns of E^{mu}_{nu} (NOT symmetric) + Esq
        if (gdd_3Dout)     ndat = ndat + 7 ! add (10-3)x independent components of g_{mu,nu} (-3 because shift=0)
        if (fluid_constraints_3Dout) ndat = ndat + 4 ! add 4 x fluid constraints - mag and escale for Mom and Ham fluid
        allocate(datas(ndat,nx,nx,nx),descriptors(ndat))

        datas(1,:,:,:)   = vel0(:,:,:,nt);   descriptors(1)  = "vel[0]"
        datas(2,:,:,:)   = vel1(:,:,:,nt);   descriptors(2)  = "vel[1]"
        datas(3,:,:,:)   = vel2(:,:,:,nt);   descriptors(3)  = "vel[2]"
        datas(4,:,:,:)   = lorentz;          descriptors(4)  = "Lorentz"
        datas(5,:,:,:)   = detg;             descriptors(5)  = "detg"
        datas(6,:,:,:)   = tracek(:,:,:,nt); descriptors(6)  = "trK"
        datas(7,:,:,:)   = Hamraw;           descriptors(7)  = "H raw"
        datas(8,:,:,:)   = sqrt(Hamescale2); descriptors(8)  = "|H| e_scale"
        datas(9,:,:,:)   = sqrt(magMom2);    descriptors(9)  = "|M| raw"
        datas(10,:,:,:)  = sqrt(Momescale2); descriptors(10) = "|M| e_scale"
        datas(11,:,:,:)  = tracer;           descriptors(11) = "traceR (slice)"
        datas(12,:,:,:)  = tracefourr;       descriptors(12) = "4R"
        if (writetilde) then
           ! Write tilde quantities; include in description
           datas(13,:,:,:)  = rho;           descriptors(13)  = "rho (tilde)"
           datas(14,:,:,:)  = sigma2;        descriptors(14)  = "sigma2 (fluid restframe; tilde)"
           datas(15,:,:,:)  = tracerfluid;   descriptors(15)  = "traceR (fluid restframe; tilde)"
           if (tderivs) then
               datas(16,:,:,:)  = theta;         descriptors(16)  = "Theta (fluid restframe; tilde)"
           else
               datas(16,:,:,:)  = theta;         descriptors(16)  = "theta (fluid restframe; tilde; linear)"
           endif
           datas(17,:,:,:)  = w2;            descriptors(17)  = "w2 (fluid restframe; tilde)"
        else
           ! Write the regular quantities
           datas(13,:,:,:)  = rho;           descriptors(13)  = "rho"
           datas(14,:,:,:)  = sigma2;        descriptors(14)  = "sigma2 (fluid restframe)"
           datas(15,:,:,:)  = tracerfluid;   descriptors(15)  = "traceR (fluid restframe)"
           if (tderivs) then
               datas(16,:,:,:)  = theta;         descriptors(16)  = "Theta (fluid restframe)"
           else
               datas(16,:,:,:)  = theta;         descriptors(16)  = "theta (fluid restframe; linear)"
           endif
           datas(17,:,:,:)  = w2;            descriptors(17)  = "w2 (fluid restframe)"
        endif
        datas(18,:,:,:)  = rhoAij;           descriptors(18) = "rho_Aij"
        datas(19,:,:,:)  = maga;             descriptors(19) = "|a|"

        ! -----------------------------------------------------------
        !
        ! ADDITIONAL OUTPUT BELOW, CANNOT HAVE ALL OF THESE AT ONCE.
        !        ALL START AT INDEX 19
        !
        ! -----------------------------------------------------------
        ni = 20 ! always the case after we have written default data
        if (shear_3Dout) then
            if (shearUd_3Dout) then
                ! We want to write sigma^{mu}_{nu} -- all 16x components
                datas(ni,:,:,:)   = sigmaUd(1,1,:,:,:);    descriptors(ni) = "sigma^0_0 (fluid restframe)"
                datas(ni+1,:,:,:) = sigmaUd(1,2,:,:,:);  descriptors(ni+1) = "sigma^0_x (fluid restframe)"
                datas(ni+2,:,:,:) = sigmaUd(1,3,:,:,:);  descriptors(ni+2) = "sigma^0_y (fluid restframe)"
                datas(ni+3,:,:,:) = sigmaUd(1,4,:,:,:);  descriptors(ni+3) = "sigma^0_z (fluid restframe)"
                datas(ni+4,:,:,:) = sigmaUd(2,1,:,:,:);  descriptors(ni+4) = "sigma^x_0 (fluid restframe)"
                datas(ni+5,:,:,:) = sigmaUd(2,2,:,:,:);  descriptors(ni+5) = "sigma^x_x (fluid restframe)"
                datas(ni+6,:,:,:) = sigmaUd(2,3,:,:,:);  descriptors(ni+6) = "sigma^x_y (fluid restframe)"
                datas(ni+7,:,:,:) = sigmaUd(2,4,:,:,:);  descriptors(ni+7) = "sigma^x_z (fluid restframe)"
                datas(ni+8,:,:,:) = sigmaUd(3,1,:,:,:);  descriptors(ni+8) = "sigma^y_0 (fluid restframe)"
                datas(ni+9,:,:,:) = sigmaUd(3,2,:,:,:);  descriptors(ni+9) = "sigma^y_x (fluid restframe)"
                datas(ni+10,:,:,:) = sigmaUd(3,3,:,:,:);  descriptors(ni+10) = "sigma^y_y (fluid restframe)"
                datas(ni+11,:,:,:) = sigmaUd(3,4,:,:,:);  descriptors(ni+11) = "sigma^y_z (fluid restframe)"
                datas(ni+12,:,:,:) = sigmaUd(4,1,:,:,:);  descriptors(ni+12) = "sigma^z_0 (fluid restframe)"
                datas(ni+13,:,:,:) = sigmaUd(4,2,:,:,:);  descriptors(ni+13) = "sigma^z_x (fluid restframe)"
                datas(ni+14,:,:,:) = sigmaUd(4,3,:,:,:);  descriptors(ni+14) = "sigma^z_y (fluid restframe)"
                datas(ni+15,:,:,:) = sigmaUd(4,4,:,:,:);  descriptors(ni+15) = "sigma^z_z (fluid restframe)"
                ! now if we have written these, ni will increase by + 16
                ni = ni + 16
            else
                ! We want to write sigma_{mu,nu} -- 10x indep. components, symmetric
                datas(ni,:,:,:)   = sigmadd(1,1,:,:,:);    descriptors(ni) = "sigma_00 (fluid restframe)"
                datas(ni+1,:,:,:) = sigmadd(1,2,:,:,:);  descriptors(ni+1) = "sigma_0x (fluid restframe)"
                datas(ni+2,:,:,:) = sigmadd(1,3,:,:,:);  descriptors(ni+2) = "sigma_0y (fluid restframe)"
                datas(ni+3,:,:,:) = sigmadd(1,4,:,:,:);  descriptors(ni+3) = "sigma_0z (fluid restframe)"
                datas(ni+4,:,:,:) = sigmadd(2,2,:,:,:);  descriptors(ni+4) = "sigma_xx (fluid restframe)"
                datas(ni+5,:,:,:) = sigmadd(2,3,:,:,:);  descriptors(ni+5) = "sigma_xy (fluid restframe)"
                datas(ni+6,:,:,:) = sigmadd(2,4,:,:,:);  descriptors(ni+6) = "sigma_xz (fluid restframe)"
                datas(ni+7,:,:,:) = sigmadd(3,3,:,:,:);  descriptors(ni+7) = "sigma_yy (fluid restframe)"
                datas(ni+8,:,:,:) = sigmadd(3,4,:,:,:);  descriptors(ni+8) = "sigma_yz (fluid restframe)"
                datas(ni+9,:,:,:) = sigmadd(4,4,:,:,:);  descriptors(ni+9) = "sigma_zz (fluid restframe)"
                ! now if we have written these, ni will increase by + 10
                ni = ni + 10
            endif
        endif
        if (Rdd_3Dout) then
            ! write R_mu,nu
            datas(ni,:,:,:)   = fourRdd(1,1,:,:,:);    descriptors(ni) = "4R_00"
            datas(ni+1,:,:,:) = fourRdd(1,2,:,:,:);  descriptors(ni+1) = "4R_0x"
            datas(ni+2,:,:,:) = fourRdd(1,3,:,:,:);  descriptors(ni+2) = "4R_0y"
            datas(ni+3,:,:,:) = fourRdd(1,4,:,:,:);  descriptors(ni+3) = "4R_0z"
            datas(ni+4,:,:,:) = fourRdd(2,2,:,:,:);  descriptors(ni+4) = "4R_xx"
            datas(ni+5,:,:,:) = fourRdd(2,3,:,:,:);  descriptors(ni+5) = "4R_xy"
            datas(ni+6,:,:,:) = fourRdd(2,4,:,:,:);  descriptors(ni+6) = "4R_xz"
            datas(ni+7,:,:,:) = fourRdd(3,3,:,:,:);  descriptors(ni+7) = "4R_yy"
            datas(ni+8,:,:,:) = fourRdd(3,4,:,:,:);  descriptors(ni+8) = "4R_yz"
            datas(ni+9,:,:,:) = fourRdd(4,4,:,:,:);  descriptors(ni+9) = "4R_zz"
            ! now if we have written these, ni will increase by + 10
            ni = ni + 10
        endif
        if (Aij_3Dout) then
            ! Write A_ij
            datas(ni,:,:,:)   = Aijdd(1,1,:,:,:);    descriptors(ni) = "A_xx"
            datas(ni+1,:,:,:) = Aijdd(1,2,:,:,:);  descriptors(ni+1) = "A_xy"
            datas(ni+2,:,:,:) = Aijdd(1,3,:,:,:);  descriptors(ni+2) = "A_xz"
            datas(ni+3,:,:,:) = Aijdd(2,2,:,:,:);  descriptors(ni+3) = "A_yy"
            datas(ni+4,:,:,:) = Aijdd(2,3,:,:,:);  descriptors(ni+4) = "A_yz"
            datas(ni+5,:,:,:) = Aijdd(3,3,:,:,:);  descriptors(ni+5) = "A_zz"
            ! now if we have written these, ni will increase by + 6
            ni = ni + 6 ! no more output; but keep this just in case
        endif
        if (EUd_3Dout) then
            ! Write E^{mu}_{nu} (Electric Weyl) -- NOt symmetric, output 16x components
            datas(ni,:,:,:)   = ElecUd(1,1,:,:,:);    descriptors(ni) = "E^0_0"
            datas(ni+1,:,:,:) = ElecUd(1,2,:,:,:);  descriptors(ni+1) = "E^0_x"
            datas(ni+2,:,:,:) = ElecUd(1,3,:,:,:);  descriptors(ni+2) = "E^0_y"
            datas(ni+3,:,:,:) = ElecUd(1,4,:,:,:);  descriptors(ni+3) = "E^0_z"
            datas(ni+4,:,:,:) = ElecUd(2,1,:,:,:);  descriptors(ni+4) = "E^x_0"
            datas(ni+5,:,:,:) = ElecUd(2,2,:,:,:);  descriptors(ni+5) = "E^x_x"
            datas(ni+6,:,:,:) = ElecUd(2,3,:,:,:);  descriptors(ni+6) = "E^x_y"
            datas(ni+7,:,:,:) = ElecUd(2,4,:,:,:);  descriptors(ni+7) = "E^x_z"
            datas(ni+8,:,:,:) = ElecUd(3,1,:,:,:);  descriptors(ni+8) = "E^y_0"
            datas(ni+9,:,:,:) = ElecUd(3,2,:,:,:);  descriptors(ni+9) = "E^y_x"
            datas(ni+10,:,:,:) = ElecUd(3,3,:,:,:);  descriptors(ni+10) = "E^y_y"
            datas(ni+11,:,:,:) = ElecUd(3,4,:,:,:);  descriptors(ni+11) = "E^y_z"
            datas(ni+12,:,:,:) = ElecUd(4,1,:,:,:);  descriptors(ni+12) = "E^z_0"
            datas(ni+13,:,:,:) = ElecUd(4,2,:,:,:);  descriptors(ni+13) = "E^z_x"
            datas(ni+14,:,:,:) = ElecUd(4,3,:,:,:);  descriptors(ni+14) = "E^z_y"
            datas(ni+15,:,:,:) = ElecUd(4,4,:,:,:);  descriptors(ni+15) = "E^z_z"
            datas(ni+16,:,:,:) = Esq;                descriptors(ni+16) = "E^2"
            ! now if we have written these, ni will increase by + 17
            ni = ni + 17
        endif
        if (gdd_3Dout) then
            ! Write g_{mu,nu} as well
            datas(ni,:,:,:)   = -alp**2;            descriptors(ni) = "g_00"
            datas(ni+1,:,:,:) = gij(1,:,:,:,nt);  descriptors(ni+1) = "g_xx"
            datas(ni+2,:,:,:) = gij(2,:,:,:,nt);  descriptors(ni+2) = "g_xy"
            datas(ni+3,:,:,:) = gij(3,:,:,:,nt);  descriptors(ni+3) = "g_xz"
            datas(ni+4,:,:,:) = gij(4,:,:,:,nt);  descriptors(ni+4) = "g_yy"
            datas(ni+5,:,:,:) = gij(5,:,:,:,nt);  descriptors(ni+5) = "g_yz"
            datas(ni+6,:,:,:) = gij(6,:,:,:,nt);  descriptors(ni+6) = "g_zz"
            ni = ni + 7
        endif
        if (fluid_constraints_3Dout) then
            ! write magMom and momescale2 in fluid frame
            datas(ni,:,:,:)   = sqrt(magFluidMom2);       descriptors(ni) = "|M| raw fluid frame"
            datas(ni+1,:,:,:) = sqrt(FluidMomescale2);  descriptors(ni+1) = "|M| e_scale fluid frame"
            datas(ni+2,:,:,:) = FluidHam;               descriptors(ni+2) = "H raw fluid frame"
            datas(ni+3,:,:,:) = sqrt(FluidHamescale2);  descriptors(ni+3) = "|H| e_scale fluid frame"
            ni = ni + 4
        endif

        call write_hdf5(nx,it,time,ndat,datas,descriptors)
        deallocate(datas,descriptors)
    endif

    !
    ! Release memory for all 3D data arrays
    !
    deallocate(tracerfluid,lorentz,sigma2,w2,theta,detg,magMom2,Momescale2,Hamraw,&
         & Hamescale2,gamijk,tracefourr,maga)
    if (shear_3Dout) then
        if (shearUd_3Dout) then
            deallocate(sigmaUd)
        else
            deallocate(sigmadd)
        endif
    endif
    if (Rdd_3Dout)   deallocate(fourRdd)
    if (Aij_3Dout)   deallocate(Aijdd)
    if (EUd_3Dout)   deallocate(ElecUd,Esq)
    if (fluid_constraints_3Dout) deallocate(magFluidMom2,FluidMomescale2,FluidHam,FluidHamescale2)

    ! --------------------------------------------------------------------------------
    !
    ! ----- WRITE AVERAGED DATA -----
    !
    ! --------------------------------------------------------------------------------

    call print_info(" Writing scalar (averaged) data ... ",loc)
    if (writeomegas) then
       call write_avg(omegam,"omegam",it,time,nspheres,rad,domain_type,&
        & "\eta, \Omega_{m}",randorigins)
       call write_avg(omegaR,"omegaR",it,time,nspheres,rad,domain_type,&
        & "\eta, \Omega_{R}",randorigins)
       call write_avg(omegaQ,"omegaQ",it,time,nspheres,rad,domain_type,&
        & "\eta, \Omega_{Q}",randorigins)
    endif

    if (writeaD) call write_avg(aDb,"aDb",it,time,nspheres,rad,domain_type,&
        & "\eta, a_{D}^{b}",randorigins)
    if (writehubble) then
        ! Write Hubble parameter:
        call write_avg(hub,"hubble",it,time,nspheres,rad,domain_type,&
            & "\eta, H_{D}",randorigins)
        if (writetilde .eqv. .False.) then
            ! write global Hubble from NON tilde theta as well
            !    --> normalise: in loop it's just <Theta>, so need to norm by \Gam and 3
            hub_allbox_notilde = hub_allbox_notilde / (3._dp * avglorentz_allbox)
            call write_avg((/ hub_allbox_notilde /),"hubble_notilde",it,time,1,&
                & tmprad,"all","\eta, H_{D} (all, no tilde)",randorigins)
        endif
    endif
    if (writedelta .and. rad/=0.) then
       call write_avg(delta,"delta",it,time,nspheres,rad,domain_type,&
        & "\eta, <rho>^b_D/bar{rho} - 1 (\delta)",randorigins)
    endif

  end subroutine compute_ricci



end module ricci
