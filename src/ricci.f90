!-------------------------------------------------------------------------
!
!  Module to compute the ricci tensor (and other things) from the metric
!
!     This module also implements the averaging scheme for general spacetime
!     foliations from Buchert, Mourier, & Roy (2019)
!
!-------------------------------------------------------------------------
module ricci
  use options, only:dp,c_double,tinitial,dit,dtfac,writeomegas,writeaD,writehubble,&
       writedelta,constraints_l1norm,tderivs,fourth,nord_dt,&
       write3D,clen,writetilde
  use roots, only:get_ricci_component,get_backreaction_omegas,fluid_restframe_curvature,&
       get_expansion_shear_vort,get_christoffel,get_fluid_curvature
  use spines, only:calc_average_within_radius,calc_average_within_radius_manyspheres
  use gauge, only:get_dtalp
  use manipulations, only:get_metric_at_pos,trace,inv3x3,get_metric_at_stencil
  use tripreports, only:write_avg,write_max,write_hdf5
  use violation, only:constraint_violation,get_L1_constraint_violation
  use init, only:initialise,get_t0_data
  use prints, only:print_info,print_error
  implicit none

contains
  !
  ! This is the main bones routine that calls other routines to calculate:
  !      - Ricci tensor and scalar (4-D and 3-D)
  !      - Christoffel symbols
  !      - trK, expansion, shear, vorticity, acceleration
  !      - perform spatial averaging, cosmological parameters and backreaction
  !      - constraint violation
  !      - write 3D output files
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

    real(c_double) :: lorentzijk,tracekijk,tracerfluidijk,tracerijk,xi,xj,xk,time
    real(c_double) :: thetaijk,sigma2ijk,sigddijk(4,4),detgijk,rhoijk,w2ijk
    real(c_double) :: xvals(nx),dt,avgrho_allbox,avgrho_allbox_notilde,avglorentz_allbox
    real(c_double) :: dtalp,alpijk,sigUUijk(4,4)
    real(c_double) :: fac,fac2,thetaavg,traceravg,rhoavg,theta2avg,sigma2avg,vort2avg,tmprad
    real(c_double), dimension(nspheres)    :: hub,Qd,omegaQ
    real(c_double), dimension(nspheres)    :: omegaR,omegam,vDbt0,aDb,rho0,delta
    real(c_double), dimension(3,3)         :: gdown,gup
    real(c_double), dimension(13,3,3)      :: gdown_atstencil,kdown_atstencil,kud_atstencil
    real(c_double), dimension(3,nspheres)  :: randorigins
    real(c_double) :: fvelU_atstencil(13,3),dialp(3),dtudt,dtudi(3),diudt(3)
    real(c_double) :: rxx,rxy,rxz,ryy,ryz,rzz,rtensor(6),threeRij(3,3)
    real(c_double) :: fluidRijk,fourRicci,fourRddijk(4,4)
    real(c_double) :: aUijk(4),magaijk,Rvsqijk,adijk(4),umudijk(4)

    integer, parameter :: navgs = 7
    real(c_double), dimension(navgs,nspheres)  :: avgs
    real(c_double), dimension(:,:,:), allocatable :: detg,theta,sigma2,w2,vsq,maga
    real(c_double), dimension(:,:,:), allocatable :: tracerfluid,tracer
    real(c_double), dimension(:,:,:), allocatable :: magMom2,Momescale2,Hamraw,Hamescale2,tracefourr
    real(c_double), dimension(:,:,:,:),     allocatable :: datas,tracek,kij_tn,gij_tn,acceld,fveld
    real(c_double), dimension(:,:,:,:,:),   allocatable :: sigmadd,fourRdd
    real(c_double), dimension(:,:,:,:,:,:), allocatable :: gamijk

    character(len=clen) :: rho0file,volt0file
    character(len=clen), allocatable :: descriptors(:)
    character(len=100) :: message,loc

    integer :: k,j,i,it,unitr0,unitv0,l,m,n,ndat
    logical :: notearly  ! if True, then we calc \Theta. if False, it is too early (i.e. we can't take time derivs yet)
    loc      = "Ricci"   ! routine location for info/errors
    notearly = .False.

    ! -----------------------------------------------------
    !
    ! Allocate memory for HDF5 output and for things we need (e.g. trK, Christoffels)
    !
    ! -----------------------------------------------------
    allocate(tracerfluid(nx,nx,nx),vsq(nx,nx,nx),sigma2(nx,nx,nx),&
         & theta(nx,nx,nx),detg(nx,nx,nx),magMom2(nx,nx,nx),Momescale2(nx,nx,nx),&
         & Hamraw(nx,nx,nx),Hamescale2(nx,nx,nx),tracer(nx,nx,nx),&
         & gamijk(3,3,3,nx,nx,nx),tracek(nx,nx,nx,nt),tracefourr(nx,nx,nx))
    allocate(kij_tn(6,nx,nx,nx),gij_tn(6,nx,nx,nx)) ! K_{ij}, g_{ij} at the current time

    ! -----------------------------------------------------
    !
    ! Intialise some grid, sphere, things etc. and filenames
    !
    ! -----------------------------------------------------
    time = times(nt) ! the current time
    call initialise(nx,nspheres,rad,time,xmin,xmax,dx,gotrho,xvals,randorigins,dt,&
         & it,notearly,volt0file,rho0file)

    ! -----------------------------------------------------
    !
    ! get initial data vol(t=0) and rho(t=0) (do we even need rho0?)
    !
    ! -----------------------------------------------------
    vDbt0 = 0._dp
    if (time/=tinitial) call get_t0_data(nspheres,rad,xmin,xmax,volt0file,rho0file,vDbt0,rho0)

    !
    ! initialise our averages
    avgrho_allbox = 0._dp;  avgrho_allbox_notilde = 0._dp; avglorentz_allbox = 0._dp
    avgs = 0._dp; omegam = 0._dp; omegaR = 0._dp; omegaQ = 0._dp
    ! A workaround for intel compilers for passing hard-coded zero radius for avgs that are always all-box
    tmprad = 0._dp

    ! -----------------------------------------------------
    !
    ! Print some user information
    !
    ! -----------------------------------------------------
    if (fourth) then
       call print_info("Using 4th order spatial derivatives ",loc)
    else
       call print_info("Using 2nd order spatial derivatives ",loc)
    endif
    if (tderivs) then
       if (notearly) then
          write(message,"(a,i1)") "Taking time derivatives of order: ",nord_dt
          call print_info(message,loc)
       else
          call print_info("Too early to take backward time derivatives. Using an approx. &
          & for expansion & shear correct to O(v^2) (for now) ",loc)
       endif
    else
       call print_info("Not taking time derivatives. Using an approx. for expansion &
       & & shear correct to O(v^2). ",loc)
    endif
    call print_info("Using avergaing formalism of Buchert+(2019) ",loc)
    !
    ! --------------------------------------------------------------------------------------------

    !
    ! Loop to get some things we need to get before the main loop
    !       -- we need trK here because we take D_i K for momentum constraint
    !       -- may as well get all Christoffels here to save MANY calls in main loop
    !
    call print_info("  Calculating Christoffel symbols and trace(K_ij) ... ",loc)
    kij_tn = kij(:,:,:,:,nt); gij_tn = gij(:,:,:,:,nt)
    !$omp parallel do default (none) &
    !$omp shared(nx,nt,dx,gij,kij,gamijk,tracek,gij_tn) &
    !$omp private(i,j,k,l,m,n,gdown)
    do k=1,nx
       do j=1,nx
          do i=1,nx
             !
             ! Loop over components of Christoffel symbols
             do l=1,3
                do m=1,3
                   do n=1,3
                      call get_christoffel(i,j,k,gij_tn,nx,dx,l,m,n,gamijk(l,m,n,i,j,k))
                   enddo
                enddo
             enddo
             !
             ! Calculate traceK
             !   -- loop over time and get whole time stencil (for 4R calc which has dt K)
             do n=1,nt
                call get_metric_at_pos(i,j,k,nx,gij(:,:,:,:,n),gdown)
                tracek(i,j,k,n) = trace(kij(:,i,j,k,n),gdown)
             enddo

          enddo
       enddo
    enddo
    !$omp end parallel do

    call print_info("  Calculating many things and averaging ... ",loc)
    !$omp parallel do default (none) &
    !$omp shared(nx,nt,dx,time,times,it,tracerfluid,rho,dens,tracek,kij,kij_tn,gij_tn) &
    !$omp shared(notearly,gotrho,gotalp,vel0,vel1,vel2,vsq,xvals) &
    !$omp shared(sigma2,w2,dt,gamijk,magMom2,Momescale2,Hamraw,Hamescale2,tmprad) &
    !$omp shared(nspheres,rad,randorigins,detg,alp,theta,tracer,sigmadd,fourRdd) &
    !$omp shared(maga,tracefourr,acceld,fveld) &
    !$omp private(i,j,k,gdown,gup,dtalp,gdown_atstencil,kdown_atstencil) &
    !$omp private(rxx,rxy,rxz,ryy,ryz,rzz,kud_atstencil,fourRicci,fourRddijk) &
    !$omp private(aUijk,adijk,magaijk,Rvsqijk,umudijk,fluidRijk) &
    !$omp private(fvelU_atstencil,dialp,dtudt,dtudi,diudt,rtensor,threeRij) &
    !$omp private(detgijk,thetaijk,sigma2ijk,sigddijk,sigUUijk,w2ijk,lorentzijk) &
    !$omp private(alpijk,fac,fac2,tracekijk,tracerfluidijk,tracerijk,rhoijk) &
    !$omp private(thetaavg,traceravg,rhoavg,theta2avg,sigma2avg,xi,xj,xk,vort2avg) &
    !$omp reduction(+:avgs,avgrho_allbox,avglorentz_allbox,vDbt0)
    do k=1,nx
       do j=1,nx
          do i=1,nx
             !
             ! reduce need to access memory throughout loop
             xi = xvals(i); xj = xvals(j); xk = xvals(k)
             !
             ! Get metric and K_ij at *spatial* finite differencing stencil
             !    and invert metric to get g^{ij}
             !
             call get_metric_at_stencil(i,j,k,nx,gij_tn,gdown_atstencil)
             gdown = gdown_atstencil(1,:,:)
             call inv3x3(gdown,gup,detgijk)
             call get_metric_at_stencil(i,j,k,nx,kij_tn,kdown_atstencil)

             !
             ! Get the relevant scalars to average
             !
             ! 1. get lapse time derivatives
             ! 2. get fluid expansion \Theta and shear \sigma^2
             !       (also w2, sigma_{mu,nu}, w_{mu,nu})
             !      -- note 2. is correct to O(v^2) if notearly=False or tderivs=False in options
             !
             alpijk    = alp(i,j,k)
             tracekijk = tracek(i,j,k,nt)
             call get_dtalp(alpijk,tracekijk,dtalp)
             call get_expansion_shear_vort(i,j,k,nx,nt,dx,dit*dt,notearly,gdown_atstencil,kdown_atstencil(1,:,:),&
                  & gamijk(:,:,:,i,j,k),alp,dtalp,gup,vel0,vel1,vel2,tracekijk,fvelU_atstencil,umudijk,dialp,&
                  & dtudt,dtudi,diudt,sigma2ijk,thetaijk,w2ijk,lorentzijk,sigddijk,sigUUijk)
             !
             ! calculate volume of each nsphere spheres at t=0
             !     -- this will be read from file for t/=tinit for a_D calc
             if (time==tinitial) call calc_average_within_radius_manyspheres(xi,xj,xk,rad,nspheres,randorigins,&
                  & detgijk,dx,1,vDbt0,(/ lorentzijk /))

             if (gotrho) then
                !
                ! we have primrho data - store this position's data
                rhoijk = rho(i,j,k)
             else
                !
                ! we have no primrho data: calculate it using conserved dens (D - see GRHydro doc)
                rhoijk     = dens(i,j,k) / ( sqrt(detgijk) * lorentzijk )
                rho(i,j,k) = rhoijk
             endif
             !
             !     Get curvature scalar in fluid rest frame \mathcal{R}
             !
             ! ----------------------------------------------------
             ! First, get 3-Ricci curvature + put into (6) for trace and (3,3) for fluidR call
             call get_ricci_component(i,j,k,gamijk,nx,dx,1,1,rxx)
             call get_ricci_component(i,j,k,gamijk,nx,dx,1,2,rxy)
             call get_ricci_component(i,j,k,gamijk,nx,dx,1,3,rxz)
             call get_ricci_component(i,j,k,gamijk,nx,dx,2,2,ryy)
             call get_ricci_component(i,j,k,gamijk,nx,dx,2,3,ryz)
             call get_ricci_component(i,j,k,gamijk,nx,dx,3,3,rzz)
             rtensor   = (/ rxx, rxy, rxz, ryy, ryz, rzz /)
             tracerijk = trace(rtensor,gdown)
             threeRij(1,1) = rxx; threeRij(2,2) = ryy; threeRij(3,3) = rzz
             threeRij(1,2) = rxy; threeRij(2,1) = rxy
             threeRij(1,3) = rxz; threeRij(3,1) = rxz
             threeRij(2,3) = ryz; threeRij(3,2) = ryz
             !
             if (tderivs .and. notearly) then
                 ! Calc the full fluidR def. in (4.15) of Buchert+(2019)
                 call get_fluid_curvature(nx,nt,dx,dit*dt,i,j,k,gdown_atstencil,gup,kij,kdown_atstencil,&
                     & tracek,alp,threeRij,gamijk(:,:,:,i,j,k),thetaijk,lorentzijk,fvelU_atstencil,dtudt,dtudi,&
                     & dtalp,dialp,diudt,kud_atstencil,aUijk,adijk,magaijk,fluidRijk,fourRicci,fourRddijk,Rvsqijk)
             else
                ! Set some things we can't calculate without time derivs to zero
                fourRicci = 0._dp; aUijk = 0._dp; adijk = 0._dp; magaijk = 0._dp
                ! Use the Hamiltonian-like constraint (4.16) insteadfor fluidR
                fluidRijk = fluid_restframe_curvature(rhoijk,thetaijk,sigma2ijk,w2ijk)
             endif
             !
             ! Build scalars to average according to new Buchert+ (2019) formalism
             !   (note: these are all multiplied by \gamma to convert to <\psi>^b_D from <\psi>_D
             !          see doc/ for more details)
             !
             fac       = alpijk / lorentzijk
             fac2      = fac**2
             thetaavg  = lorentzijk * fac * thetaijk
             traceravg = lorentzijk * fac2 * fluidRijk
             rhoavg    = lorentzijk * fac2 * rhoijk
             sigma2avg = lorentzijk * fac2 * sigma2ijk
             theta2avg = lorentzijk * fac2 * thetaijk**2
             vort2avg  = lorentzijk * fac2 * w2ijk

             !
             ! Calculate averages within spheres or over whole domain
             !   avgs: theta, R, rho, sigma2, theta2, w2, lorentz
             !
             call calc_average_within_radius_manyspheres(xi,xj,xk,rad,nspheres,randorigins,&
                  & detgijk,dx,navgs,avgs,(/ thetaavg, traceravg, rhoavg, sigma2avg, theta2avg, &
                  & vort2avg, lorentzijk /))

             !
             ! Calculate average density over the WHOLE BOX; for density contrast \delta
             !     and V^b_all = \sum_all \Gamma for global avg density, delta
             !  NOTE: all box avg ==> rad=0., nspheres=1, randorigins=0.
             !
             call calc_average_within_radius(xi,xj,xk,(/tmprad,tmprad,tmprad/),tmprad,detgijk,dx,&
                  & avgrho_allbox,rhoavg)
             call calc_average_within_radius(xi,xj,xk,(/tmprad,tmprad,tmprad/),tmprad,detgijk,dx,&
                  & avglorentz_allbox,lorentzijk)
             !
             ! Calculate constraint violation
             !
             call constraint_violation(i,j,k,nx,dx,gdown_atstencil,kdown_atstencil,gamijk,&
                  & tracek(:,:,:,nt),rhoijk,vel0(i,j,k,nt),vel1(i,j,k,nt),vel2(i,j,k,nt),lorentzijk,&
                  & tracerijk,magMom2(i,j,k),Momescale2(i,j,k),Hamraw(i,j,k),Hamescale2(i,j,k))
             !
             ! Store some things for 3D output
             !
             vsq(i,j,k)         = 1._dp - 1._dp / lorentzijk**2 ! this is \gamma_{ij} v^i v^j
             detg(i,j,k)        = detgijk
             tracefourr(i,j,k)  = fourRicci
             tracer(i,j,k)      = tracerijk
             if (writetilde) then
                ! Write tilde 3D quantities
                rho(i,j,k)         = fac2 * rhoijk ! overwrite this value with tilde rho
                theta(i,j,k)       = fac * thetaijk
                sigma2(i,j,k)      = fac2 * sigma2ijk
                w2(i,j,k)          = fac2 * w2ijk
                tracerfluid(i,j,k) = fac2 * fluidRijk
             else
                ! Write regular quantities
                ! do nothing for rho since we have this stored already
                theta(i,j,k)       = thetaijk
                sigma2(i,j,k)      = sigma2ijk
                w2(i,j,k)          = w2ijk
                tracerfluid(i,j,k) = fluidRijk
             endif

          enddo
       enddo
    enddo
    !$omp end parallel do

    ! --------------------------------------------------------------------------------
    !   Get backreaction terms and cosmological params
    !
    !      NOTE: Volume of comoving domain is:
    !      V^b_D = V_D * <lorentz>_D
    !            = \int_D lorentz \sqrt{gamma} d^3X (1/V_D cancels -- see "user guide")
    !
    call print_info("  Calculating backreaction terms and cosmological parameters ... ",loc)
    call get_backreaction_omegas(it,time,rad,nspheres,randorigins,avgs,navgs,avgrho_allbox,&
         & avgrho_allbox_notilde,avglorentz_allbox,vDbt0,Qd,hub,omegam,omegaR,omegaQ,delta,aDb)
    if (nspheres<10) then
        ! Print some of these if the output won't be too much
        write(message,"(a,ES14.7)") "      <rho>_D(all) = ",avgrho_allbox
        call print_info(message,loc)
        write(message,"(a,ES14.7)") "          V_D(all) = ",avglorentz_allbox
        call print_info(message,loc)
    endif

    ! --------------------------------------------------------------------------------
    ! Get L1/L2 error for constraints -- and writes to files (scalars)
    ! --------------------------------------------------------------------------------
    if (constraints_l1norm) call get_L1_constraint_violation(nx,it,time,Hamraw,Hamescale2,magMom2,Momescale2)

    ! --------------------------------------------------------------------------------
    ! Write t_initial rho, volume
    ! --------------------------------------------------------------------------------
    if (time==tinitial) then
       call print_info(" Writing t = t_init data: "//trim(rho0file)//', '//trim(volt0file),loc)
       open(newunit=unitr0,file=rho0file,status='replace')
       open(newunit=unitv0,file=volt0file,status='replace')
       rho0  = avgs(3,:) ! avg rho
       write(unitr0,*) rho0
       write(unitv0,*) vDbt0
       close(unitr0)
       close(unitv0)
    endif

    ! --------------------------------------------------------------------------------
    !
    !     Write 3D data to HDF5 files
    !
    ! --------------------------------------------------------------------------------
    if (write3D) then
       call print_info(" Writing 3D (HDF5) data ... ",loc)
       ! number of 3-D arrays we want to write
       ndat = 17
       allocate(datas(ndat,nx,nx,nx),descriptors(ndat))

       datas(1,:,:,:)   = vel0(:,:,:,nt);   descriptors(1)  = "vel[0]"
       datas(2,:,:,:)   = vel1(:,:,:,nt);   descriptors(2)  = "vel[1]"
       datas(3,:,:,:)   = vel2(:,:,:,nt);   descriptors(3)  = "vel[2]"
       datas(4,:,:,:)   = vsq;              descriptors(4)  = "vsq"
       datas(5,:,:,:)   = detg;             descriptors(5)  = "detg"
       datas(6,:,:,:)   = tracek(:,:,:,nt); descriptors(6)  = "trK"
       datas(7,:,:,:)   = Hamraw;           descriptors(7)  = "H raw"
       datas(8,:,:,:)   = sqrt(Hamescale2); descriptors(8)  = "|H| e_scale"
       datas(9,:,:,:)   = sqrt(magMom2);    descriptors(9)  = "|M| raw"
       datas(10,:,:,:)  = sqrt(Momescale2); descriptors(10) = "|M| e_scale"
       datas(11,:,:,:)  = tracer;           descriptors(11) = "traceR (slice)"
       datas(12,:,:,:)  = tracefourr;       descriptors(12) = "4R"
       if (writetilde) then
          ! Write tilde quantities; include appropriate description
          datas(13,:,:,:)  = rho;           descriptors(13)  = "rho (tilde)"
          datas(14,:,:,:)  = sigma2;        descriptors(14)  = "sigma2 (fluid restframe; tilde)"
          datas(15,:,:,:)  = tracerfluid;   descriptors(15)  = "traceR (fluid restframe; tilde)"
          datas(16,:,:,:)  = theta;         descriptors(16)  = "theta (fluid restframe; tilde)"
          datas(17,:,:,:)  = w2;            descriptors(17)  = "w2 (fluid restframe; tilde)"
       else
          ! Write the regular quantities
          datas(13,:,:,:)  = rho;           descriptors(13)  = "rho"
          datas(14,:,:,:)  = sigma2;        descriptors(14)  = "sigma2 (fluid restframe)"
          datas(15,:,:,:)  = tracerfluid;   descriptors(15)  = "traceR (fluid restframe)"
          datas(16,:,:,:)  = theta;         descriptors(16)  = "theta (fluid restframe)"
          datas(17,:,:,:)  = w2;            descriptors(17)  = "w2 (fluid restframe)"
       endif
       !
       ! Write the data & free memory
       call write_hdf5(nx,it,time,ndat,datas,descriptors)
       deallocate(datas,descriptors)
    endif
    !
    ! Free memory for all 3D data arrays
    deallocate(tracerfluid,vsq,sigma2,w2,theta,detg,magMom2,Momescale2,Hamraw,&
         & Hamescale2,gamijk,tracer,tracefourr,tracek,kij_tn,gij_tn)

    ! --------------------------------------------------------------------------------
    !
    !   Write averaged data
    !
    ! --------------------------------------------------------------------------------

    call print_info(" Writing scalar (averaged) data ... ",loc)
    if (writeomegas) then
       call write_avg(omegam,"omegam",it,time,nspheres,rad,"\eta, \Omega_{m}",randorigins)
       call write_avg(omegaR,"omegaR",it,time,nspheres,rad,"\eta, \Omega_{R}",randorigins)
       call write_avg(omegaQ,"omegaQ",it,time,nspheres,rad,"\eta, \Omega_{Q}",randorigins)
    endif
    if (writeaD)     call write_avg(aDb,"aDb",it,time,nspheres,rad,"\eta, a_{D}^{b}",randorigins)
    if (writehubble) call write_avg(hub,"hubble",it,time,nspheres,rad,"\eta, H_{D}",randorigins)
    if (writedelta .and. rad/=0.) then
       call write_avg(delta,"delta",it,time,nspheres,rad,"\eta, <rho>^b_D/bar{rho} - 1 (\delta)",randorigins)
    endif

  end subroutine compute_ricci

end module ricci
