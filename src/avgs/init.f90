module init
  !
  ! A module to initialise some things for Mescaline
  !
  !        - General init routine
  !        - Make filenames for spheres and (old) 3D ascii data
  !        - Read in t = 0 data for volume, density
  !
  use options, only:pi,dtfac,c_double,tinitial,dit,&
      tderivs,sphere_rseed,looprad,radmax,nord_dt,clen,domain_type,&
      read_domain_origins,origins_fname
  use manipulations, only:trace
  use prints, only:print_error,print_info
  use tripreports, only:filename_avg
  use random, only:init_random_seed
  use analytic_solns, only:get_init_rho
  implicit none

contains

  !
  ! Initialise some things:
  !      - grid
  !      - read prev. timestep files (if required)
  !      - setup and define spheres for averagine (if required)
  !
  subroutine initialise(nx,nspheres,rad,time,xmin,xmax,dx,gotrho,xvals,randorigins,&
       dt,it,notearly)
    integer, intent(in) :: nx,nspheres
    real(c_double), intent(in) :: rad,xmin,xmax,dx,time
    logical, intent(in) :: gotrho

    real(c_double), intent(out) :: xvals(nx),dt,randorigins(3,nspheres)
    integer, intent(out) :: it
    logical, intent(inout) :: notearly

    real(c_double) :: maxorigin,minorigin,randnums(3,nspheres),t1,t2,t3,t4
    real(c_double) :: t_thresh
    integer :: l,i,ounit

    character(len=100) :: message,originsfile
    character(len=6) :: loc
    loc = "Init" ! the location to be printed to screen for this routine

    ! ---------------------------------------------------
    !
    ! Some general User info
    !
    ! ---------------------------------------------------
    call print_info(" Initialising...",loc)

    write(message,"(a,f5.3)") "    t_init = ",tinitial
    call print_info(message,loc)

    if (tderivs) then
       write(message,"(a,i3,a)")'    Skipping every ',dit,' iterations '
       call print_info(message,loc)
    endif

    if (gotrho .eqv. .False.) call print_info("     No primrho data: will calculate in each loop ",loc)

    ! ----------------------------------------------------
    !
    ! Set x-value array, timestep dt and current iteration from time
    !
    ! ----------------------------------------------------

    do l=1,nx
       xvals(l) = xmin + (l-1)*dx
    enddo
    dt = dtfac * dx
    it = nint((time-tinitial)/dt) ! this can be just *below* the it we want when dt is small, so need nint()

    ! ----------------------------------------------------
    !
    ! Calculate current iteration and previous times we need (if applicable)
    !      --> and find out if we are at a late enough time to take time derivatives
    !          depending on nord_dt as set in options
    !
    ! ----------------------------------------------------
    !
    ! t1 is the first time AFTER the initial conditions that we have output
    !     we need this to calculate \Theta -- specifically a time derivative of u_t
    !
    !     --> NOTE: dit is the spacing between HDF5 output files NOT simulation timestep ***
    !
    t_thresh = 9999.
    if (tderivs) then
         if (nord_dt<2) call print_error(" Minimum required nord_dt = 2. Please re-set and re-run. ",2,loc)
         t1 = tinitial + dit * dt         ! t_i + dt: minumum req. is two previous timesteps; for 2nd order, always need these
         t2 = tinitial + 2.d0 * dit * dt ! t_i + 2*dt
         ! set threshold time after which we can take time derivatives
         t_thresh = t2

         if (nord_dt >= 3) then
            !
            ! We want at least a 3rd order scheme, so need minimum 4 steps total (including t_n)
            !
            t3 = tinitial + 3.d0 * dit * dt ! t_i + 3*dt
            ! set threshold time after which we can take time derivatives
            t_thresh = t3

            if (nord_dt >= 4) then
               !
               ! We want at least a 4th order scheme, so need minumum 5 steps total (incluidng t_n)
               !
               t4 = tinitial + 4.d0 * dit * dt ! t_i + 4 * dt
               ! set threshold time after which we can take time derivatives
               t_thresh = t4

               if (nord_dt > 4) call print_error(" nord_dt > 4 NOT IMPLEMENTED. Please set nord_dt <= 4 and re-run. ",2,loc)
            endif
        endif
        notearly = .False.
        if (time >= t_thresh) notearly = .True.

    endif

    ! ----------------------------------------------------
    !
    ! Build random origins for averaging // or read from a file
    !
    ! ----------------------------------------------------

    if (rad==0.) then
        !
        ! We want the whole box average; no subdomains
        !
        call print_info("    Averaging over the WHOLE BOX ",loc)
        randorigins = 0.d0
    elseif (rad>0. .and. nspheres==1) then
        !
        ! Set origin to be centre for single sphere -- can be changed in future
        !
        randorigins = (xmax + xmin) / 2.d0
        write(message,"(a,f10.5)") "    Averaging over ONE "//trim(domain_type)//" with radius = ",rad !," Mpc "
        call print_info(message,loc)
        write(message,"(a,f10.5,a,f10.5,a,f10.5)") "    With origin: x = ",randorigins(1,1),", y = ",randorigins(2,1),&
                ", z = ",randorigins(3,1)
        call print_info(message,loc)

    else
        !
        ! We want many spheres; first check if we want to read them from file or not
        !
        if (looprad) then
            maxorigin = xmax - radmax  ! All sphere origins should be in the same range
            minorigin = xmin + radmax
        else
            maxorigin = xmax - rad  ! if randnum = 1 --> we want this number
            minorigin = xmin + rad  ! if randnum = 0 --> we want this number
        endif
        ! Give some info
        write(message,"(a,i6,a,f10.5)") "    Averaging over ",nspheres," "//trim(domain_type)//"s with radii = ",rad !," Mpc "
        call print_info(message,loc)

        if (read_domain_origins) then
            !
            ! Read the origins from a file
            !
            call print_info("    With many different origins read from: "//trim(origins_fname),loc)
            ! open the file
            open(file=origins_fname,newunit=ounit,status='old')
            do i=1,nspheres
                read(ounit,*) randorigins(:,i)
            enddo
            close(ounit)
            ! check if any of the origins are too close to the boundary
            if (any(randorigins>maxorigin) .or. any(randorigins<minorigin)) then
                call print_info(" Some of your chosen spheres will overlap the boundary (I cant deal with that yet). Fix it pls.",loc)
            endif
        else
            !
            ! Draw the random numbers for the sphere origins
            !
            call init_random_seed(sphere_rseed)
            call RANDOM_NUMBER(randnums)
            randorigins = (1.d0 - randnums) * minorigin + randnums * maxorigin
            !
            ! Write origins of spheres to file -- if we didnt read from one
            write(originsfile,"(a,i4.4,a,i4.4,a)") trim(domain_type)//'_origins_r',int(rad),'_nsph',nspheres,'.dat'
            open(file=originsfile,newunit=ounit,status='replace')
            write(ounit,*) "# x, y, z (radius = ",rad,", nspheres = ",nspheres,")"
            do i=1,nspheres
                write(ounit,*) randorigins(:,i)
            enddo
            close(ounit)
            call print_info("    With many different origins written to: "//trim(originsfile),loc)
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
  !  ** deprecated **  // leaving here for funsies
  !
  subroutine get_t0_data(nspheres,rad,xmin,xmax,f1,f2,data1,data2)
    integer, intent(in) :: nspheres
    real(c_double), intent(in) :: rad,xmin,xmax
    character(len=clen), intent(in) :: f1, f2

    real(c_double), intent(out) :: data1(nspheres),data2(nspheres)! volt0, rho0

    real(c_double) :: boxlen,rhoinitial,hubinit,rhostar
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
         select case(domain_type)
         case("sphere")
             ! Set an approximate spherical volume
             ! this is precise to ~0.0006 for highest res sims (256,1G and 128,500)
             data1 = (4.d0 / 3.d0) * pi * rad**3
         case("cube")
             ! Set an approximate cubic volume
             ! (I haven't tested how precise this is but it should be similar to above)
             data1 = rad**3
         case default
             call print_error("Please choose a valid domain_type",2,loc)
         end select

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
       call get_init_rho(rhoinitial,hubinit,rhostar)
       data2 = rhoinitial
    endif

  end subroutine get_t0_data



end module init
