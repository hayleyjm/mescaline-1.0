!---------------------------------------------------------------------
!
! Mescaline: Interesting things extracted from cactus
! Code to read in node-separated ET HDF5 files and output trip reports
!  No analysis done here; modified version of main code
! Written by Hayley Macpherson &  Daniel Price, 2017-2020
!
!---------------------------------------------------------------------
program mescaline
  use cactoidae, only:read_cactus_file
  use prints,    only:print_mescaline,print_info,print_error,print_compliment
  use options,   only:c_double,clen,tinitial
  use tripreports, only:write_hdf5
  use tools_rw, only:get_it

 implicit none
 integer :: i,nargs,nx,it,ndat
 real(c_double) :: dx,time,xmin(3),xmax(3)
 !
 ! arrays to store the data
 real(c_double), allocatable, dimension(:,:,:,:) :: kij,gij
 real(c_double), allocatable, dimension(:,:,:)   :: alp,vel0,vel1,vel2,rho,dens
 real(c_double), allocatable, dimension(:,:,:,:) :: datas

 character(len=clen), allocatable :: descriptors(:)
 character(len=clen) :: message,loc,filename
 loc = "Mescaline"

 call print_mescaline()
 print*, "  Reading ET files and re-writing output into trip reports. Doing no analysis. "

 !
 ! get number of filenames from command line
 !
 nargs = command_argument_count()

 if (nargs >= 1) then
    write(message,"(a,i3,a)") "  --> re-writing a total of ",nargs," files. "

    do i=1,nargs

      ! -------------------------------------------------------------------
      !
      ! We do NOT want to calculate time derivatives; use one timestep ALWAYS
      !
      ! -------------------------------------------------------------------
      call get_command_argument(i,filename)

      ! Read this file and get the grid data for this time
      call read_cactus_file(1,filename,time,dx,xmin,xmax,nx,gij,&
      & kij,alp,vel0,vel1,vel2,rho,dens)
      it = get_it(dx,time)
      print*, ' it = ',it,' time = ',time,' tinit = ',tinitial

      !
      ! Now we want to write the data to the trip reports
      !
      ndat = 17
      allocate(datas(ndat,nx,nx,nx),descriptors(ndat))

      datas(1,:,:,:)   = vel0;   descriptors(1)  = "vel[0]"
      datas(2,:,:,:)   = vel1;   descriptors(2)  = "vel[1]"
      datas(3,:,:,:)   = vel2;   descriptors(3)  = "vel[2]"

      datas(4,:,:,:) = gij(1,:,:,:);  descriptors(4) = "g_xx"
      datas(5,:,:,:) = gij(2,:,:,:);  descriptors(5) = "g_xy"
      datas(6,:,:,:) = gij(3,:,:,:);  descriptors(6) = "g_xz"
      datas(7,:,:,:) = gij(4,:,:,:);  descriptors(7) = "g_yy"
      datas(8,:,:,:) = gij(5,:,:,:);  descriptors(8) = "g_yz"
      datas(9,:,:,:) = gij(6,:,:,:);  descriptors(9) = "g_zz"

      datas(10,:,:,:) = kij(1,:,:,:);  descriptors(10) = "K_xx"
      datas(11,:,:,:) = kij(2,:,:,:);  descriptors(11) = "K_xy"
      datas(12,:,:,:) = kij(3,:,:,:);  descriptors(12) = "K_xz"
      datas(13,:,:,:) = kij(4,:,:,:);  descriptors(13) = "K_yy"
      datas(14,:,:,:) = kij(5,:,:,:);  descriptors(14) = "K_yz"
      datas(15,:,:,:) = kij(6,:,:,:);  descriptors(15) = "K_zz"

      datas(16,:,:,:)   = alp;   descriptors(16)  = "alpha"
      datas(17,:,:,:)   = rho;   descriptors(17)  = "rho"

      call write_hdf5(nx,it,time,ndat,datas,descriptors)

      ! free up the space
      deallocate(gij,kij,vel0,vel1,vel2,rho,dens)
      deallocate(datas,descriptors)
      write(*,"(a)") ""
    enddo

 else
    print "(/,a,/)",' Usage: mescaline *.hdf5 '
    stop
 endif

 write(*,"(a)") ""
 call print_info(" All done. You made it. See the trip report/s for details. ")
 call print_compliment()

end program mescaline
