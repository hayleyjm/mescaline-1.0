!-------------------------------------------------------------------------
!
!  Module that handles the front end to the data read
!  and rearranges the data into a sensible grid
!
!-------------------------------------------------------------------------
module cactoidae
 use, intrinsic :: iso_c_binding, only:c_int
 use options, only:c_double,clen
 use prints,  only:print_info,print_error
 implicit none

 integer, parameter :: ignoretl = 1
 integer :: ncells
 logical :: gotalp,gotrho ! gotrho(gotalp) true if we have primitive rho (lapse) data location in HDF5 file

 real(c_double), allocatable :: dat(:,:)
 integer,        allocatable :: iamtype(:)

contains


  !-------------------------------------------------------------------------
  !
  !  Main routine to read cactus data file and pass data off for analysis
  !
  !-------------------------------------------------------------------------
  subroutine read_cactus_file(n,filename,time,dx,xmin,xmax,nx,gij,kij,alp,vel0,vel1,vel2,rho,dens)
    use cactushdf5read, only:open_cactus_hdf5_file,close_cactus_hdf5_file,blocklabel,read_cactus_hdf5_data
    use asciiutils,     only:cstring
    use prints,         only:print_info,print_error

    integer, intent(in) :: n ! where we are in the list of files being run (purely for printing)
    character(len=*), intent(in) :: filename

    real(c_double), intent(out) :: dx,time,xmin(3),xmax(3)
    real(c_double), allocatable, dimension(:,:,:,:), intent(out) :: gij,kij ! (6,nx,nx,nx)
    real(c_double), allocatable, dimension(:,:,:),   intent(out) :: vel0,vel1,vel2,rho,dens,alp ! (nx,nx,nx)
    integer, intent(out) :: nx

    integer :: ncells,ncol,nsteps,ndim,ndimV,istep,ierr
    character(len=120) :: message,loc
    loc = "Cactoid" ! current location for message printing

    write(message,"(a,i1)") 'Reading '//trim(filename)//' ignore tl = ',ignoretl
    call print_info(message,loc)
    !
    ! Read header
    !
    istep = 1
    call open_cactus_hdf5_file(cstring(filename),istep,ncells,ncol,nsteps,ndim,ndimV,time,ignoretl,ierr)
    !if (debug) print*, ' nf = ',nf,' file = ',filename,'ncells = ',ncells,'ncol = ',ncol,'time = ',time,'ierr = ',ierr
    if (ierr /= 0) then
       write(message,"(a)") ' ERROR: could not read from '//trim(filename)
       call close_cactus_hdf5_file(ierr)
       call print_error(message,2,loc)
       return
    endif

    !
    ! Allocate memory for dat array (for this timestep)
    !
    write(message,"(a,i10,a,i5,a,i3,a,f10.3,f10.3,f10.3,f10.3)") ' >>> Allocating: ncells = ',ncells,', ncol = ',&
         ncol,', nsteps = ',nsteps,', time = ',time
    call print_info(message,loc)
    allocate(dat(ncells,ncol))
    if (allocated(iamtype) .eqv. .False.) allocate(iamtype(ncells))
    dat = 0.

    blocklabel(1:5) = (/'x ','y ','z ','dx','m '/)
    !
    ! Read data from file
    !
    istep = 1
    call read_cactus_hdf5_data(cstring(filename),istep,ncells,time,dx,ignoretl,ierr)
    ! after the call above, dat is filled with the hdf5 data

    !
    ! extract metric, kij, fluid vars from dat array
    !
    ! these data arrays are allocated INSIDE extract_data. deallocate them after we have stored them
    call extract_data(ncells,ncol,n,iamtype,dx,xmin,xmax,blocklabel,nx,gij,&
      & kij,alp,vel0,vel1,vel2,rho,dens)

    deallocate(dat)
    call close_cactus_hdf5_file(ierr)

  end subroutine read_cactus_file




  !
  ! this module extracts the metric and extrinsic curvature
  ! from the dataset and regrids it so we can compute derivatives
  !
  subroutine extract_data(ncells,ncol,n,itype,dx,xmin,xmax,labels,nx,gij,kij,&
        & alp,vel0,vel1,vel2,rho,dens)
    use cactushdf5read, only:find_metric
    integer,          intent(in) :: ncells,ncol,n ! 'n' is where we are in time loop in mescaline.f90
    real(c_double),   intent(in) :: dx
    integer,          intent(in) :: itype(ncells)
    character(len=*), intent(in) :: labels(ncol)

    real(c_double), allocatable, intent(out) :: gij(:,:,:,:),kij(:,:,:,:)
    real(c_double), allocatable, dimension(:,:,:), intent(out) :: alp,vel0,vel1,vel2,rho,dens
    real(c_double), intent(out) :: xmin(3),xmax(3)
    integer, intent(out) :: nx

    integer :: i,ii,jj,kk
    integer :: igxx,igxy,igxz,igyy,igyz,igzz
    integer :: ikxx,ikxy,ikxz,ikyy,ikyz,ikzz
    integer :: irho,ialp,ivel0,ivel1,ivel2
    real(c_double) :: xi,yi,zi
    character(len=clen) :: message,loc
    loc = " extract_data" ! our current location for print messages

    ! count number of actual cells (not ghost cells)
    nx = nint(count(itype==1)**(1./3.))

    ! find location of metric in columns
    call find_metric(ncol,labels,igxx,igxy,igxz,igyy,igyz,igzz,&
         ikxx,ikxy,ikxz,ikyy,ikyz,ikzz,irho,ialp,&
         ivel0,ivel1,ivel2,gotrho,gotalp)
    !
    ! allocate memory for individual data
    !
    if (allocated(alp)) deallocate(alp) ! safeguard because we don't deallocate this every time in cactoid
    if (allocated(gij)) deallocate(gij)
    allocate(gij(6,nx,nx,nx),alp(nx,nx,nx),kij(6,nx,nx,nx),vel0(nx,nx,nx),vel1(nx,nx,nx),&
         vel2(nx,nx,nx),rho(nx,nx,nx),dens(nx,nx,nx))

    ! find grid size
    do i=1,3
       xmin(i) = minval(dat(:,i)) + 3.*dx
       xmax(i) = maxval(dat(:,i)) - 2.*dx
    enddo
    if (n==1) then
      ! Only want to print grid data once, since it's always the same!
      write(message,"(a,i8)") " >>>  Grid data:    nx = ",nx
      call print_info(message,loc)
      write(message,"(a,f11.5)") " >>>                dx = ",dx
      call print_info(message,loc)
      write(message,"(a,f12.1,f12.1,f12.1)") " >>>              xmin = ",xmin
      call print_info(message,loc)
      write(message,"(a,f12.1,f12.1,f12.1)") " >>>              xmax = ",xmax
      call print_info(message,loc)
    endif
    write(message,"(a)") " Regridding and extracting metric... "
    call print_info(message,loc)

    do i=1,ncells
       if (itype(i)==1) then ! not for ghost cells

          xi = dat(i,1)
          yi = dat(i,2)
          zi = dat(i,3)
          ii = get_cell_index(xi,xmin(1),dx,nx)
          jj = get_cell_index(yi,xmin(2),dx,nx)
          kk = get_cell_index(zi,xmin(3),dx,nx)

          gij(1,ii,jj,kk) = dat(i,igxx)
          gij(2,ii,jj,kk) = dat(i,igxy)
          gij(3,ii,jj,kk) = dat(i,igxz)
          gij(4,ii,jj,kk) = dat(i,igyy)
          gij(5,ii,jj,kk) = dat(i,igyz)
          gij(6,ii,jj,kk) = dat(i,igzz)
          if (gotalp) alp(ii,jj,kk) = dat(i,ialp)

          kij(1,ii,jj,kk) = dat(i,ikxx)
          kij(2,ii,jj,kk) = dat(i,ikxy)
          kij(3,ii,jj,kk) = dat(i,ikxz)
          kij(4,ii,jj,kk) = dat(i,ikyy)
          kij(5,ii,jj,kk) = dat(i,ikyz)
          kij(6,ii,jj,kk) = dat(i,ikzz)
          vel0(ii,jj,kk)  = dat(i,ivel0)
          vel1(ii,jj,kk)  = dat(i,ivel1)
          vel2(ii,jj,kk)  = dat(i,ivel2)
          if (gotrho) then
             rho(ii,jj,kk)  = dat(i,irho)
          else
             dens(ii,jj,kk) = dat(i,irho)
          endif
       endif
    enddo

  end subroutine extract_data



  !
  ! find cell index corresponding to a given position
  !
  integer function get_cell_index(x,xmin,dx,nx) result(i)
    real(c_double), intent(in) :: x,xmin,dx
    integer,        intent(in) :: nx
    character(len=clen) :: loc,message
    loc = " get_cell_index"

    i = nint((x - xmin)/dx) + 1
    if (i < 1) then
      write(message,"(a,i3)") " i < 1. i = ",i
      call print_error(message,2,loc)
    endif
    if (i > nx) then
      write(message,"(a,i3)") " i > nx. i = ",i
      call print_error(message,2,loc)
    endif

  end function get_cell_index



!-------------------------------------------------------------------------
!
!  The following routines are callback routines called by the c
!  utilities
!
!-------------------------------------------------------------------------
subroutine read_cactus_hdf5_data_fromc(icol,ntot,np,temparr) bind(c)
  integer(kind=c_int), intent(in) :: icol,ntot,np
  real(kind=c_double), intent(in) :: temparr(np)
  integer(kind=c_int) :: icolput
  integer :: nmax,i1,i2

  icolput = icol
  i1 = ntot-np+1
  i2 = ntot

  !if (debugmode) print "(a,i2,a,i8,a,i8)",&
  !'DEBUG: reading column ',icol,' -> '//trim(label(icolput))//' parts ',i1,' to ',i2

  ! check column is within array limits
  if (icolput.gt.size(dat(1,:)) .or. icolput.eq.0) then
     print "(a,i2,a)",' ERROR: column = ',icolput,' out of range in receive_data_fromc'
     return
  endif
  if (i2 > size(dat(:,1))) then
     print*,' ERROR with index range: ',i1,':',i2,' exceeds size ',size(dat(:,1)),' for column ',icol
     read*
     return
  endif
  ! ensure no array overflows
  nmax = min(i2,size(dat(:,1)))

  ! copy data into main array
  dat(i1:i2,icolput) = real(temparr(1:np))

  return
end subroutine read_cactus_hdf5_data_fromc

subroutine read_cactus_itype_fromc(ntot,np,itype) bind(c)
  integer(kind=c_int), intent(in) :: ntot,np
  integer(kind=c_int), intent(in) :: itype(np)
  integer :: i1,i2,len_type

  i1 = ntot-np+1
  i2 = ntot
  ! set particle type
  len_type = size(iamtype(:))
  if (len_type.gt.1) then
     if (i2 > len_type) then
        print*,'error with itype length',i2,len_type
        return
     endif
     iamtype(i1:i2) = itype(1:np)
  endif

  return
end subroutine read_cactus_itype_fromc

end module cactoidae
