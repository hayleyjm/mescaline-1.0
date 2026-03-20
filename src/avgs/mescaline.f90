!---------------------------------------------------------------------
!
! Mescaline: Interesting things extracted from cactus
! Code to perform post-simulation analysis on Cactus HDF5 data files
! Written by Hayley Macpherson &  Daniel Price, 2017-2020
!
!---------------------------------------------------------------------
program mescaline
  use cactoidae, only:read_cactus_file
  use analysis,  only:send_for_analysis
  use prints,    only:print_mescaline,print_info,print_error,print_compliment
  use options,   only:tderivs,nord_dt,restart,c_double,clen
 implicit none
 integer :: i,nargs,n,idx,nfiles,nfiles_toread
 integer :: nx
 real(c_double) :: dx,xmin(3),xmax(3)
 real(c_double), allocatable :: times(:)
 !
 ! arrays to store the data to send for analysis - at a particular time
 real(c_double), allocatable, dimension(:,:,:,:) :: kijn,gijn
 real(c_double), allocatable, dimension(:,:,:)   :: alp,vel0n,vel1n,vel2n,rho,dens
 ! - at all times
 real(c_double), allocatable, dimension(:,:,:,:,:) :: kij,gij
 real(c_double), allocatable, dimension(:,:,:,:) :: vel0,vel1,vel2!,rho,dens

 character(len=100), allocatable :: filenames(:)
 character(len=clen) :: message,loc
 logical :: first,keepdata
 loc = "Mescaline"

 call print_mescaline()

 !
 ! Number of files we need if taking time derivs
 !   --> NOTE: this is distinct from the number of files we will ACTUALLY read
 !              (which is situation-dependent, see below)
 !
 nfiles_toread = nord_dt + 1

 !
 ! get number of filenames from command line
 !
 nargs = command_argument_count()

 !
 ! Safeguard for in case we're running with only one file
 !   as argument but want time derivs:
 !
 if (nargs==1 .and. tderivs .eqv. .True.) then
    call print_error("WARNING: only one file in argument list, cannot take time derivatives. &
        & Please adjust your run script OR set tderivs = False.",2,loc)
 endif

 first    = .True.    ! flag for the first time we have i>=nfiles_toread
 keepdata = .False.   ! whether to store previous 4 timesteps in dat_keep
 !
 ! a dummy allocate to remove some compile-time warnings...
 allocate(gij(1,1,1,1,1),kij(1,1,1,1,1),vel0(1,1,1,1),vel1(1,1,1,1),&
    & vel2(1,1,1,1))!,rho(1,1,1,1),dens(1,1,1,1))

 if (nargs >= 1) then

    do i=1,nargs

       if (tderivs) then
          ! -------------------------------------------------------------------
          !
          ! We want to calculate time derivatives; we need multiple timesteps
          !
          ! -------------------------------------------------------------------
          if (i >= nfiles_toread) then
              !
              ! We now have enough files to grab (nord_dt-1) previous files & this one
              !
              keepdata = .True. ! We want to keep previous timesteps for next time from now on!

              nfiles = nfiles_toread

              write(message,"(a,i2,a)") "Reading ",nfiles," timesteps ... "
              call print_info(message,loc)
              allocate(filenames(nfiles))
              !
              ! Read all filenames from command line (should be no problem doing this every time?)
              !
              do n=1,nfiles
                   idx = i - (nfiles - n)
                   call get_command_argument(idx,filenames(n))
              enddo
              !
              ! Loop over files and read (all if first, only most recent if not)
              if (first) then
                  !
                  first = .False. ! only do this once
                  !
                  ! re-allocate times (may have prev. been only one time step)
                  if (allocated(times)) deallocate(times)
                  allocate(times(nfiles))
                  !
                  ! read ALL nfiles files in
                  !
                  do n=1,nfiles
                      ! Read this file and get the grid data for this time
                      call read_cactus_file(n,filenames(n),times(n),dx,xmin,xmax,&
                        & nx,gijn,kijn,alp,vel0n,vel1n,vel2n,rho,dens)
                      !
                      ! Store this data in multi-time arrays
                      if (n==1) then
                        ! allocate memory for all-time arrays ONCE
                        !
                        ! check for allocation first, since if we run first few files individually
                        !      then these will have already been allocated...
                        if (allocated(kij)) deallocate(kij,gij,vel0,vel1,vel2)!,rho,dens)
                        allocate(gij(6,nx,nx,nx,nfiles),kij(6,nx,nx,nx,nfiles),vel0(nx,nx,nx,nfiles),&
                             & vel1(nx,nx,nx,nfiles),vel2(nx,nx,nx,nfiles))!,rho(nx,nx,nx,nfiles),dens(nx,nx,nx,nfiles))
                     endif
                     gij(:,:,:,:,n) = gijn
                     kij(:,:,:,:,n) = kijn
                     vel0(:,:,:,n)  = vel0n
                     vel1(:,:,:,n)  = vel1n ! I HATE this. We do it 3 times. Weird to put it into a subroutine
                     vel2(:,:,:,n)  = vel2n !    but maybe there's some restructure we can do to so we can avoid this!?
                     !rho(:,:,:,n)   = rhon
                     !dens(:,:,:,n)  = densn
                      !
                      ! deallocate memory for temporary n arrays (others deallocated in extract_data)
                      deallocate(gijn,kijn,vel0n,vel1n,vel2n)!,rhon,densn)
                  enddo

              else ! NOT First

                  !
                  ! We've already read all other timesteps in, only read in most recent time
                  !
                  call read_cactus_file(nfiles,filenames(nfiles),times(nfiles),dx,xmin,&
                    & xmax,nx,gijn,kijn,alp,vel0n,vel1n,vel2n,rho,dens)
                  !
                  ! We have already shifted the previous times' data to keep it
                  gij(:,:,:,:,nfiles) = gijn
                  kij(:,:,:,:,nfiles) = kijn
                  vel0(:,:,:,nfiles)  = vel0n
                  vel1(:,:,:,nfiles)  = vel1n
                  vel2(:,:,:,nfiles)  = vel2n
                  !rho(:,:,:,nfiles)   = rhon
                  !dens(:,:,:,nfiles)  = densn

                  ! deallocate memory for temporary n arrays (others deallocated in extract_data)
                  deallocate(gijn,kijn,vel0n,vel1n,vel2n)!,rhon,densn)
              endif

          else ! i < nfiles_toread

              nfiles = 1
              if (restart) then

                   !
                   ! This run is a restart of some previous run:
                   !    --> user will have passed in nord_dt+1 previous files so we can
                   !        continue time derivatives as normal.
                   !    --> so, we already have output for the first nfiles files, and we
                   !        can ignore them.
                   !
                   if (i==1) then
                     write(message,"(a,i2,a)") " RESTART: ignoring first ",nord_dt," files "
                     call print_info(message,loc)
                   endif
                   cycle

              else ! NOT a restart run

                   !
                   ! We need to run a few single-file runs to get things going
                   !
                   if (i==1) then
                     write(message,"(a,i2,a)") "Reading in one timestep until I can read ",nfiles_toread," timesteps "
                     call print_info(message,loc)
                   endif
                   if (allocated(times)) deallocate(times) ! a good check to have - not sure if necessary here?
                   allocate(filenames(nfiles),times(nfiles))
                   call get_command_argument(i,filenames(1))
                   !
                   ! Read this file and get the grid data for this time
                   call read_cactus_file(1,filenames(1),times(1),dx,xmin,xmax,nx,&
                    & gijn,kijn,alp,vel0n,vel1n,vel2n,rho,dens)
                   !
                   ! Store this data in multi-time arrays

                   ! allocate memory for all-time arrays
                   !
                   ! if we come into this loop more than once, need to deallocate
                   if (allocated(kij)) deallocate(gij,kij,vel0,vel1,vel2)!,rho,dens)
                   allocate(gij(6,nx,nx,nx,nfiles),kij(6,nx,nx,nx,nfiles),vel0(nx,nx,nx,nfiles),&
                      & vel1(nx,nx,nx,nfiles),vel2(nx,nx,nx,nfiles))!,rho(nx,nx,nx,nfiles),dens(nx,nx,nx,nfiles))

                   gij(:,:,:,:,1) = gijn
                   kij(:,:,:,:,1) = kijn
                   vel0(:,:,:,1)  = vel0n
                   vel1(:,:,:,1)  = vel1n
                   vel2(:,:,:,1)  = vel2n
                   !rho(:,:,:,1)   = rhon
                   !dens(:,:,:,1)  = densn

                   ! deallocate memory for temporary n arrays (others deallocated in extract_data)
                   deallocate(gijn,kijn,vel0n,vel1n,vel2n)!,rhon,densn)
              endif

          endif

       else
          ! -------------------------------------------------------------------
          !
          ! We do NOT want to calculate time derivatives; use one timestep ALWAYS
          !
          ! -------------------------------------------------------------------
          if (i==1) call print_info("Reading one file at a time")
          nfiles = 1
          if (allocated(times)) deallocate(times)
          allocate(filenames(nfiles),times(nfiles))
          call get_command_argument(i,filenames(1))

          ! Read this file and get the grid data for this time
          call read_cactus_file(1,filenames(1),times(1),dx,xmin,xmax,nx,gijn,&
            & kijn,alp,vel0n,vel1n,vel2n,rho,dens)
          !
          ! Store this data in multi-time arrays
          !
          ! allocate memory for all-time arrays
          if (allocated(kij)) deallocate(gij,kij,vel0,vel1,vel2)!,rho,dens)
          allocate(gij(6,nx,nx,nx,nfiles),kij(6,nx,nx,nx,nfiles),vel0(nx,nx,nx,nfiles),&
          & vel1(nx,nx,nx,nfiles),vel2(nx,nx,nx,nfiles))!,rho(nx,nx,nx,nfiles),dens(nx,nx,nx,nfiles))

          gij(:,:,:,:,1) = gijn
          kij(:,:,:,:,1) = kijn
          vel0(:,:,:,1)  = vel0n
          vel1(:,:,:,1)  = vel1n
          vel2(:,:,:,1)  = vel2n
          !rho(:,:,:,1)   = rhon
          !dens(:,:,:,1)  = densn

          ! deallocate memory for temporary n arrays (others deallocated in extract_data)
          deallocate(gijn,kijn,vel0n,vel1n,vel2n)!,rhon,densn)
       endif

       ! -------------------------------------------------------------------
       !
       ! Send data in for analysis and release memory we don't need till next loop
       !
       ! -------------------------------------------------------------------
       call print_info(" Sending data for analysis ... ",loc)
       call send_for_analysis(nx,nfiles,nargs,times,xmax,xmin,dx,gij,kij,alp,vel0,vel1,vel2,rho,dens)

       !
       ! Keep some data for next loop, if we want, before deallocating
       !
       if (keepdata) then
          ! SHIFT the data backwards in time dimension rather than having a whole other array
          gij(:,:,:,:,1:nfiles-1) = gij(:,:,:,:,2:nfiles)
          kij(:,:,:,:,1:nfiles-1) = kij(:,:,:,:,2:nfiles)
          vel0(:,:,:,1:nfiles-1)  = vel0(:,:,:,2:nfiles)
          vel1(:,:,:,1:nfiles-1)  = vel1(:,:,:,2:nfiles)
          vel2(:,:,:,1:nfiles-1)  = vel2(:,:,:,2:nfiles)
          !rho(:,:,:,1:nfiles-1)   = rho(:,:,:,2:nfiles)
          !dens(:,:,:,1:nfiles-1)  = dens(:,:,:,2:nfiles)
          times(1:nfiles-1)       = times(2:nfiles)
       endif
       deallocate(filenames)
       write(*,"(a)") ""
    enddo

 else
    print "(/,a,/)",' Usage: mescaline *.hdf5 '
    stop
 endif

 deallocate(kij,vel0,vel1,vel2,rho,dens,gij,alp,times)

 write(*,"(a)") ""
 call print_info(" All done. You made it. See the trip report/s for details. ")
 call print_compliment()

end program mescaline
