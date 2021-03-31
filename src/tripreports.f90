module tripreports
  !
  ! A module to write output data for mescaline
  !   this module contains ALL output routines for mescaline
  !
  use options, only:c_double,dit,fourth,clen,tinitial
  use prints, only:print_info,print_error
  implicit none
  character(len=100) :: message
  character(len=100), parameter :: loc = " tripreports"

contains


  !
  ! A function to create filenames on-the-go for:
  !     --> 3D output data / constraints
  !      "simple" because it's just the ident and iteration
  !
  function filename_simp(ident,it)
    character(len=100) :: filename_simp
    character(len=*) :: ident
    integer :: it

    write(filename_simp,"(a,i6.6,a)") ident//'_it',it,'.dat'

  end function filename_simp



  !
  ! A function to create filenames on-the-go for:
  !     --> averages over spheres / all box
  !
  function filename_avg(ident,it,rad,nspheres)
    character(len=100) :: filename_avg
    character(len=*)   :: ident
    integer :: it,nspheres,intrad
    real(c_double) :: rad,nonintrad,rinttrad
    !
    ! Check if we are doing a full box average or spheres
    !    and name file accordingly
    !
    if (rad==0. .and. nspheres==1) then
       !
       ! We're doing a whole-box average
       !
       write(filename_avg,"(a,i6.6,a)")   ident//'_all_',it,'.dat'
    else
       !
       ! We're averaging over spheres!
       !     --> First: check if we have an integer radius, or if we need to format specially
       !
       intrad   = int(rad)        ! integer part of radius
       rinttrad = real(intrad)    ! real.. integer part of radius... for comparison
       if (rinttrad==rad) then
          ! The radius is a whole number, e.g. 100.0 or 50.0
          !    --> format with just an integer for rad
          write(filename_avg,"(a,i4.4,a,i4.4,a,i6.6,a)") ident//'_r',intrad,'_nsph',nspheres,'_it',it,'.dat'
       else
          ! There is some non-integer part of our radius, e.g. 10.5 or 6.7
          !    --> format with a real number for rad (pad with zeros in a weird way, cos Fortran...)
          nonintrad = rad - rinttrad
          write(filename_avg,"(a,i4.4,f0.2,a,i4.4,a,i6.6,a)") ident//'_r',intrad,nonintrad,'_nsph',nspheres,'_it',it,'.dat'
       endif
    endif

  end function filename_avg



  !
  ! a subroutine to write HDF5 output data
  !
  ! specifically for 3D data with dimensions (nx,nx,nx)
  !
  subroutine write_hdf5(nx,it,time,ndat,datas,descriptors,filename_amend)
    use hdf5
    integer, intent(in) :: nx,it,ndat                                  ! current iteration and number of arrays we are writing
    real(c_double), intent(in) :: time                              ! the current coordinate time
    real(c_double), intent(in) :: datas(ndat,nx,nx,nx)   ! an array of the 3D data to write to file
    character(len=*), dimension(ndat) :: descriptors                ! a char. array of the names of the data
    character(len=*), optional :: filename_amend        ! optional additional string to add to filename to e.g. distinguish from analytic vs. numerical cases

    !
    ! various things we need to open the file/dataspaces/datasets to do with the data itself
    !
    !         (types here are to correspond with the C types we need)
    !
    character(len=100) :: filename                  ! the name of the final hdf5 file
    integer(HID_T) :: file_id                      ! the ID of the HDF5 file
    integer(HID_T) :: dset_ids(ndat),dspace_id     ! the ID of the ndat datasets, and the single dataspace for all ndat

    !
    ! various things specific to the global attributes dataset
    !
    character(len=18), parameter :: attr_dsetname = "Global Attributes"  ! the name of the attr dataset
    integer, parameter :: nattrs = 4                                     ! the number of attrs we want to attach to the file
    integer(HID_T) :: attrspace_id,attrset_id,attr_ids(nattrs)           ! the IDs for the attr dataspace, dataset, and attrs themselves
    integer(HID_T) :: attrtypes(nattrs),aint_type_id,achar_type_id,arl_type_id ! the IDs for the attr types: integer, character, real
    character(len=6), dimension(nattrs) :: attrnames                     ! the names of the attributes themselves
    integer(HSIZE_T), dimension(1), parameter :: adims = (/1/)           ! the dimensions of the dataspace for the global attrs.
    integer(HSIZE_T), dimension(1), parameter :: attr_data_dims = (/1/)  ! the dimensions for the attribute data itself
    integer(SIZE_T),  parameter               :: fourthlen = 5           ! the length of the character attribute "fourth"
    character(len=fourthlen)                  :: fourthval               ! the character value of fourth = "True" or "False"

    !
    ! some things about the data itself to store in the datasets
    !
    integer, parameter :: drank = 3                                      ! rank of the datasets
    integer(HSIZE_T), dimension(drank) :: dims ! Dataset dimensions

    integer :: err,n,idx

    ! ============================================================
    !
    ! Initialize HDF5 FORTRAN interface.
    !
    ! ============================================================
    call h5open_f(err)

    ! ============================================================
    !
    ! Create a new file with default properties
    !
    ! ============================================================

    ! ------------------------------------------
    ! Write the filename based on current iteration
    ! ------------------------------------------
    if (present(filename_amend)) then
      write(filename,"(a,i6.6,a)")   "trip_report_3D_"//trim(filename_amend)//"_it",it,".h5"
    else
      write(filename,"(a,i6.6,a)")   "trip_report_3D_it",it,".h5"
    endif
    call print_info(" Writing to "//trim(filename),loc)

    ! ------------------------------------------
    ! Create the file
    ! ------------------------------------------
    !
    ! note: the flag H5F_ACC_TRUNC_F overwrites an existing file
    !          H5F_ACC_EXCL_F  to exit if file already exists
    !
    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, err)
    if (err/=0) then
       print*, ' ERROR creating file ', filename,'. STOPPING '
       stop
    endif


    ! ============================================================
    !
    ! Create dataspace, datasets and write the 3D data itself
    !
    ! ============================================================

    ! ------------------------------------------
    ! Create the dataspace
    ! ------------------------------------------
    !
    ! note all datasets here can use the same dataspace
    !   except for the attributes, which are added below
    !
    dims = (/ nx, nx, nx /)
    call h5screate_simple_f(drank, dims, dspace_id, err)

    ! ------------------------------------------
    ! Create and write to each dataset in this dataspace
    ! ------------------------------------------
    do n=1,ndat

       ! Create the dataset
       call h5dcreate_f(file_id, descriptors(n), H5T_NATIVE_DOUBLE, dspace_id, &
            dset_ids(n), err)

       ! Write to the dataset
       call h5dwrite_f(dset_ids(n), H5T_NATIVE_DOUBLE, datas(n,:,:,:), dims, err)

       ! Close the dataset
       call h5dclose_f(dset_ids(n), err)

    enddo

    ! ------------------------------------------
    ! Terminate access to the data space
    ! ------------------------------------------
    call h5sclose_f(dspace_id, err)

    ! ============================================================
    !
    ! Create the ``Global Attributes'' dataspace, dataset, and attributes
    !
    ! ============================================================


    ! ------------------------------------------
    ! Write the names of the global attributes
    ! ------------------------------------------
    !
    ! (currently can't think of a better way to do this than to hard-wire this part)
    !
    attrnames(1) = "nx"
    attrnames(2) = "time"
    attrnames(3) = "dit"
    attrnames(4) = "fourth"

    ! ------------------------------------------
    ! Make a scalar dataspace to store attrs
    ! ------------------------------------------
    !
    ! (first parameter here is: rank)
    !
    call h5screate_simple_f(1, adims, attrspace_id, err)

    ! ------------------------------------------
    ! Make an empty dataset for the global parameters
    ! ------------------------------------------
    call h5dcreate_f(file_id, attr_dsetname, H5T_NATIVE_INTEGER, attrspace_id, &
         attrset_id, err)

    ! ------------------------------------------
    ! Create datatypes for the INTEGER and CHARACTER attributes
    ! ------------------------------------------
    !
    ! (and define the length for char = true or false)
    !
    ! integers
    call h5tcopy_f(H5T_NATIVE_INTEGER, aint_type_id, err)
    ! characters
    call h5tcopy_f(H5T_NATIVE_CHARACTER, achar_type_id, err)
    call h5tset_size_f(achar_type_id, fourthlen, err)
    ! reals
    call h5tcopy_f(H5T_NATIVE_DOUBLE, arl_type_id, err)
    !
    ! ------------------------------------------
    ! Set types corresponding to attrnames above
    ! ------------------------------------------
    attrtypes(1) = aint_type_id
    attrtypes(2) = arl_type_id
    attrtypes(3) = aint_type_id
    attrtypes(4) = achar_type_id

    ! ------------------------------------------
    ! Create dataset attributes
    ! ------------------------------------------
    !
    !  note: written in reverse-order due to the way it's stored in the HDF5 file (want time at the top)
    !
    do n=1,nattrs
       idx = nattrs - (n-1)
       call h5acreate_f(attrset_id, attrnames(idx), attrtypes(idx), attrspace_id, attr_ids(idx), err)
    enddo

    ! ------------------------------------------
    ! Write the attribute data
    ! ------------------------------------------
    !
    !  (again hard-wired for now, since we can't mix types in Fortran arrays)
    !
    ! Write parameter ``fourth'' as a character rather than logical
    if (fourth) then
       fourthval = "True"
    else
       fourthval = "False"
    endif
    call h5awrite_f(attr_ids(4), achar_type_id, fourthval, attr_data_dims, err)
    call h5awrite_f(attr_ids(3), aint_type_id, (/ dit /), attr_data_dims, err)
    call h5awrite_f(attr_ids(2), arl_type_id, (/ time /), attr_data_dims, err)
    call h5awrite_f(attr_ids(1), aint_type_id, (/ nx /), attr_data_dims, err)

    ! ------------------------------------------
    ! Close the attributes
    ! ------------------------------------------
    do n=1,nattrs
       idx = nattrs - (n-1)
       call h5aclose_f(attr_ids(idx), err)
    enddo

    ! ------------------------------------------
    ! Close the attribute datatypes
    ! ------------------------------------------
    call h5tclose_f(aint_type_id, err)
    call h5tclose_f(achar_type_id, err)
    call h5tclose_f(arl_type_id, err)

    ! ------------------------------------------
    ! End access to the dataset and release resources used by it
    ! ------------------------------------------
    call h5dclose_f(attrset_id, err)

    ! ------------------------------------------
    ! Terminate access to the data space
    ! ------------------------------------------
    call h5sclose_f(attrspace_id, err)


    ! ============================================================
    !
    ! Close the file.
    !
    ! ============================================================
    call h5fclose_f(file_id, err)

    ! ============================================================
    !
    ! Close FORTRAN HDF5 interface.
    !
    ! ============================================================
    call h5close_f(err)

  end subroutine write_hdf5



  !
  ! a subroutine to write an average quantity to a file, i.e. over time
  ! * quantity can be a single number, or a list of numbers (header should reflect this)
  ! * can also pass multiple quantities to multiple files, this will assume size(quantity)=size(filename)
  ! and write them respectively, i.e. quantity(1) goes into filename(1), with header(1)
  !
  subroutine write_avg(quantity,ident,it,time,nspheres,radius,header,randorigins)
    integer, intent(in) :: nspheres,it
    real(c_double), intent(in) :: quantity(nspheres), time, radius
    character(len=*), intent(in) :: ident ! identifier of whatever we are writing
    character(len=*), intent(in), optional :: header
    real(c_double),   intent(in), optional :: randorigins(3,nspheres)

    character(len=clen) :: filename
    integer :: unit

    filename = filename_avg(ident,it,radius,nspheres)

    open(newunit=unit,file=filename,status='replace')
    call print_info(" Writing to "//trim(filename),loc)
    if (time==tinitial) then   ! We want to write randorigins and header to file
       if (radius/=0.) then    ! We have one sphere with rad/=0.
          if (present(randorigins)) then
              write(unit,"(a,i4.4,a)") '# coordinate, origin(1), origin(2), ... origin(',nspheres,')'
              write(unit,*) '# x', randorigins(1,:)
              write(unit,*) '# y', randorigins(2,:)
              write(unit,*) '# z', randorigins(3,:)
          endif
          if (present(header)) write(unit,"(a,a,a,i4.4,a)") '# ',header,'(1 ...',nspheres,')'
       else                    ! we have rad=0. (whole box avg)
          write(unit,"(a)") '# whole box average '
          if (present(header)) write(unit,"(a,a,a,i4.4,a)") '# ',header
       endif
    endif
    write(unit,*) time, quantity
    close(unit)

  end subroutine write_avg


  !
  ! a subroutine to write the MAXIMUM (abs) value of a 3D quantity to a file
  !
  subroutine write_max(nx,quantity,ident,it,time)
    integer, intent(in) :: nx,it
    real(c_double), intent(in) :: quantity(nx,nx,nx),time
    character(len=*), intent(in) :: ident ! string NAME identifying quantity to write
    real(c_double) :: maxquantity
    character(len=clen) :: filename
    integer :: unit

    !
    ! write the filename
    filename = filename_simp(trim(ident//'_max'),it)
    !
    ! calculate the MAX - of abs
    maxquantity = maxval(abs(quantity))

    !
    ! open the file and write to it
    open(newunit=unit,file=filename,status='replace')
    call print_info(" Writing to "//trim(filename),loc)
    if (it==0) then
       !
       ! Write header to file
       write(unit,"(a)") '# time, maximum |'//ident//'|'
    endif
    write(unit,*) time, maxquantity
    close(unit)

  end subroutine write_max


end module tripreports
