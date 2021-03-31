module init
  !
  ! A module to initialise some things for Mescaline
  !
  !        - General init routine
  !        - Make filenames for spheres and (old) 3D ascii data
  !        - Read in t = 0 data for volume, density
  !
  use options, only:ainit,rhoinit,pi,dtfac,c_double,tinitial,dit,dp,&
       tderivs,sphere_rseed,nord_dt,clen
  use manipulations, only:trace
  use prints, only:print_error,print_info
  use tripreports, only:filename_avg
  use random, only:init_random_seed
  implicit none

contains

  !
  ! Initialise some things:
  !      - grid
  !      - read prev. timestep files (if required)
  !      - setup and define spheres for averagine (if required)
  !
  subroutine initialise(nx,nspheres,rad,time,xmin,xmax,dx,gotrho,xvals,randorigins,&
       dt,it,notearly,volt0file,rho0file)
    integer, intent(in) :: nx,nspheres
    real(c_double), intent(in) :: rad,xmin,xmax,dx,time
    logical, intent(in) :: gotrho

    real(c_double), intent(out) :: xvals(nx),dt,randorigins(3,nspheres)
    integer, intent(out) :: it
    logical, intent(inout) :: notearly
    character(len=clen), intent(out) :: volt0file,rho0file ! filenames for some tinitial things we need

    real(c_double) :: maxorigin,minorigin,randnums(3,nspheres),t1,t2,t3,t4
    real(c_double) :: t_thresh
    integer :: l,i,ounit

    character(len=100) :: message,originsfile
    character(len=6) :: loc
    loc = "Init" ! the location to be printed to screen for this routine

    !
    ! Some general User info
    !
    call print_info(" Initialising...",loc)
    write(message,"(a,f5.3)") "    t_init = ",tinitial
    call print_info(message,loc)
    if (tderivs) then
       write(message,"(a,i3,a)")'        dit = ',dit
       call print_info(message,loc)
    endif
    if (gotrho .eqv. .False.) call print_info("     No primrho data: will calculate in each loop ",loc)

    !
    ! Set x-value array, timestep dt and current iteration from time
    !
    do l=1,nx
       xvals(l) = xmin + (l-1)*dx
    enddo
    dt = dtfac * dx
    it = nint((time-tinitial)/dt) ! this can be just *below* the it we want when dt is small, so need nint()

    !
    ! Set filenames for t=tinitial data we need if averaging
    !
    rho0file  = filename_avg('rhoavg_t0',0,rad,nspheres) ! setting it=0 here
    volt0file = filename_avg('volume_t0',0,rad,nspheres) ! since these are always for it=0 data

    !
    ! Calculate current iteration and previous times we need (if applicable)
    !      --> and find out if we are at a late enough time to take time derivatives
    !          depending on nord_dt as set in options
    !
    ! ----------------------------------------------------
    ! t1 is the first time AFTER the initial conditions that we have output
    !     we need this to calculate \Theta -- specifically a time derivative of u_t
    !     --> NOTE: dit is the spacing between HDF5 output files NOT simulation timestep ***
    !
    t_thresh = 9999.
    if (tderivs) then
       if (nord_dt<2) call print_error(" Minimum required nord_dt = 2. Please re-set and re-run. ",2,loc)
       t1 = tinitial + dit * dt         ! t_i + dt: minumum req. is two previous timesteps; for 2nd order, always need these
       t2 = tinitial + 2._dp * dit * dt ! t_i + 2*dt
       ! set threshold time after which we can take time derivatives
       t_thresh = t2

       if (nord_dt >= 3) then
          !
          ! We want at least a 3rd order scheme, so need minimum 4 steps total (including t_n)
          !
          t3 = tinitial + 3._dp * dit * dt ! t_i + 3*dt
          ! set threshold time after which we can take time derivatives
          t_thresh = t3

          if (nord_dt >= 4) then
             !
             ! We want at least a 4th order scheme, so need minumum 5 steps total (incluidng t_n)
             !
             t4 = tinitial + 4._dp * dit * dt ! t_i + 4 * dt
             ! set threshold time after which we can take time derivatives
             t_thresh = t4

             if (nord_dt > 4) call print_error(" nord_dt > 4 NOT IMPLEMENTED. Please set 2 < nord_dt<= 4 ",2,loc)
          endif
       endif
       notearly = .False.
       if (time >= t_thresh) notearly = .True.

    endif

    ! ----------------------------------------------------
    !
    ! Build random origins for averaging
    !
    ! ----------------------------------------------------
    if (rad==0.) then
       call print_info("    Averaging over the WHOLE BOX ",loc)
       randorigins = 0._dp
    elseif (rad>0. .and. nspheres==1) then
       !
       ! Set origin to be centre for single sphere
       !
       randorigins = (xmax + xmin) / 2._dp
       write(message,"(a,f7.1)") "    Averaging over ONE sphere with radius = ",rad
       call print_info(message,loc)
       write(message,"(a,f7.2,a,f7.2,a,f7.2)") "    With origin: x = ",randorigins(1,1),", y = ",randorigins(2,1),&
            ", z = ",randorigins(3,1)
       call print_info(message,loc)

    else
       !
       ! We have many spheres
       !     -- ensure spheres don't overlap boundaries
       !     -- draw random numbers as origins and re-scale to these min/max values
       !
       maxorigin = xmax - rad  ! if randnum = 1 --> we want this number
       minorigin = xmin + rad  ! if randnum = 0 --> we want this number
       call init_random_seed(sphere_rseed)
       call RANDOM_NUMBER(randnums)
       randorigins = (1._dp - randnums) * minorigin + randnums * maxorigin
       write(message,"(a,i6,a,f7.1)") "    Averaging over ",nspheres," spheres with radii = ",rad
       call print_info(message,loc)

       if (time==tinitial) then
          !
          ! Write origins of spheres to file
          !
          write(originsfile,"(a,i4.4,a,i4.4,a)") 'sphere_origins_r',int(rad),'_nsph',nspheres,'.dat'
          open(file=originsfile,newunit=ounit,status='replace')
          write(ounit,*) "# x, y, z (radius = ",rad,", nspheres = ",nspheres,")"
          do i=1,nspheres
             write(ounit,*) randorigins(:,i)
          enddo
          call print_info("    With many different origins written to: "//trim(originsfile),loc)
       else
          call print_info("    With many different origins. ",loc)
       endif
    endif

    call print_info(" Done initialising. ",loc)

  end subroutine initialise




  !
  ! Get data for volume (data1) and density (data2) at t = tinitial
  !
  !    --> Looks for a file which will be output by mescaline when run on t = tinitial data
  !           if the file doesn't exist: we display an error message, and continue with approximate values
  !
  subroutine get_t0_data(nspheres,rad,xmin,xmax,f1,f2,data1,data2)
    integer, intent(in) :: nspheres
    real(c_double), intent(in) :: rad,xmin,xmax
    character(len=clen), intent(in) :: f1, f2

    real(c_double), intent(out) :: data1(nspheres),data2(nspheres)

    real(c_double) :: boxlen
    logical :: exist1,exist2
    integer :: u1,u2,i,ierr
    character(len=150) :: message,loc
    loc = " get_t0_data"

    boxlen = (xmax - xmin)

    inquire(file=f1,exist=exist1)
    inquire(file=f2,exist=exist2)
    if (exist1) then
       !
       ! We have a (t=0) file. Read volt0 data
       !
       open(newunit=u1,file=f1,status='old',action='read',iostat=ierr)
       i=0
       do while (ierr>0)
          if (i==0) call print_info(" Waiting to open volt0 file... ",loc)
          !
          ! File exists but we cannot open it for some reason
          !     this is likely because you may be running mescaline on multiple files at once,
          !     and another run is currently reading from the file.
          !
          !     So wait in this loop until whatever program is done using it! (same for rho0)
          !
          close(u1)
          open(newunit=u1,file=f1,status='old',action='read',iostat=ierr)
          i = i+1
       enddo
       call print_info(" Reading vol(t=0) ... "//f1,loc)
       read(u1,*) data1
       close(u1)
    else
       !
       ! We DON'T have a (t=0) file: SET APPROXIMATE VOL(t=init)
       !
       write(message,"(a,f7.1,a,i6,a)") "You have not run mescaline at t = tinitial for this file using radius = ",&
            rad,", nspheres = ",nspheres,". Setting APPROXIMATE initial volume."
       call print_error(message,1,loc)
       if (rad==0.) then
          ! volume of entire box
          data1 = boxlen**3
       else
          ! volume of averaging domains
          data1 = (4._dp / 3._dp) * pi * rad**3 ! this is precise to ~0.0006 for highest res sims (256,1G and 128,500)
       endif
    endif

    if (exist2) then
       !
       ! We have a (t=0) file. Read volt0 data
       !
       open(newunit=u2,file=f2,status='old',action='read',iostat=ierr)
       i=0
       do while (ierr>0)
          if (i==0) call print_info(" Waiting to open rho0 file... ",loc)
          close(u2)
          open(newunit=u2,file=f2,status='old',action='read',iostat=ierr)
          i = i+1
       enddo
       call print_info(" Reading rho(t=0) ... "//f2,loc)
       read(u2,*) data2
       close(u2)
    else
       !
       ! We DON'T have a (t=0) file: SET rho0 as set in options.f90
       !
       write(message,"(a,f7.1,a,i6,a)") "You have not run mescaline at t = tinitial for this file using radius = ",&
            rad,", nspheres = ",nspheres,". Setting initial density as in options.f90"
       call print_error(message,1,loc)
       data2 = rhoinit
    endif

  end subroutine get_t0_data





  !
  ! Set filenames for 3D violation data -- called from constraint_violation ONLY
  !
  !    --> NOTE: this form of this routine called from violation routine that has its own spatial loop
  !                 updated version is called from INSIDE the main ricci.f90 loop
  !
  subroutine set_filenames_3D_violation_old(it,hamrelfile,momrelfile,hamrawfile,momrawfile,MomL1rawfile,&
       MomL1relfile,HamL1rawfile,HamL1relfile,h1file,h2file,M1relfile,M2relfile,M3relfile,&
       M1rawfile,M2rawfile,M3rawfile,HamrawMaxfile,HamrelMaxfile,&
       MomrawMaxfile,MomrelMaxfile)
    integer, intent(in) :: it
    character(len=*), intent(inout) :: momrelfile,momrawfile,MomL1rawfile,MomL1relfile
    character(len=*), intent(inout) :: hamrelfile,hamrawfile,h1file,h2file,HamL1rawfile,HamL1relfile
    character(len=*), intent(inout) :: M1relfile,M2relfile,M3relfile
    character(len=*), intent(inout) :: M1rawfile,M2rawfile,M3rawfile
    character(len=*), intent(inout) :: HamrawMaxfile,HamrelMaxfile,MomrawMaxfile,MomrelMaxfile

    !
    ! write raw AND rel constraints -- doesn't matter if we just name ALL files...
    !

    !
    ! 3D data filenames
    write(hamrawfile,"(a,i6.6,a)")   'Hamraw3D_it',it,'.dat'
    write(momrawfile,"(a,i6.6,a)")   'Momraw3D_it',it,'.dat'
    write(hamrelfile,"(a,i6.6,a)")   'Hamrel3D_it',it,'.dat'
    write(momrelfile,"(a,i6.6,a)")   'Momrel3D_it',it,'.dat'

    !
    ! L1 error (scalar) filenames
    write(HamL1rawfile,"(a,i6.6,a)") 'Ham_L1error_raw_it',it,'.dat'
    write(MomL1rawfile,"(a,i6.6,a)") 'Mom_L1error_raw_it',it,'.dat'
    write(HamL1relfile,"(a,i6.6,a)") 'Ham_L1error_rel_it',it,'.dat'
    write(MomL1relfile,"(a,i6.6,a)") 'Mom_L1error_rel_it',it,'.dat'

    !
    ! Max valued (scalar) filenames
    write(HamrawMaxfile,"(a,i6.6,a)") 'Ham_Max_raw_it',it,'.dat'
    write(HamrelMaxfile,"(a,i6.6,a)") 'Ham_Max_rel_it',it,'.dat'
    write(MomrawMaxfile,"(a,i6.6,a)") 'Mom_Max_raw_it',it,'.dat'
    write(MomrelMaxfile,"(a,i6.6,a)") 'Mom_Max_rel_it',it,'.dat'

    !
    ! 3D data for individual M_i
    write(M1rawfile,"(a,i6.6,a)")    'M1raw3D_it',it,'.dat'
    write(M2rawfile,"(a,i6.6,a)")    'M2raw3D_it',it,'.dat'
    write(M3rawfile,"(a,i6.6,a)")    'M3raw3D_it',it,'.dat'
    write(M1relfile,"(a,i6.6,a)")    'M1rel3D_it',it,'.dat'
    write(M2relfile,"(a,i6.6,a)")    'M2rel3D_it',it,'.dat'
    write(M3relfile,"(a,i6.6,a)")    'M3rel3D_it',it,'.dat'

    !
    ! 3D data for individual terms that make up H
    write(h1file,"(a,i6.6,a)")       'Ham_Kudud_3D_it',it,'.dat'
    write(h2file,"(a,i6.6,a)")       'Ham_K2_3D_it',it,'.dat'

  end subroutine set_filenames_3D_violation_old



end module init
