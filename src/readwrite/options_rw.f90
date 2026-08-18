module options
  !
  ! A module to set options for run-time in Mescaline (readwrite)
  !
  !  READ ALL COMMENTS CAREFULLY to
  !    ensure all parameters here are set to desired/correct values before compiling
  !
  !       (when changing any of these you must *clean* and re-compile)
  !
  use, intrinsic :: iso_c_binding, only:c_double
  implicit none

  !integer,        parameter :: dp = 8            ! precision of floats
  integer, parameter :: clen = 100               ! default length of character variables for all characters

  ! -----------------------------------------------------------------------------
  !
  ! Specify some things about the NR simulation we are analysing
  !    (most of these can be found in the .par file of the simulation)
  !
  ! ** CHANGE THESE to ensure the iteration is correct in the output files
  !
  ! -----------------------------------------------------------------------------
  real(c_double), parameter :: tinitial = 2.6168031077368945     ! initial time
  real(c_double), parameter ::    dtfac = 0.1                     ! Courant factor: dt = dtfac * dx

  ! -----------------------------------------------------------------------------
  !
  ! DO NOT change the below -- they are passed as dummy values to routines also
  !      used by other mescaline codes
  !
  ! -----------------------------------------------------------------------------
  ! Spacing (in iterations) between HDF5 files, e.g. how often we have 3D dumps
  integer, parameter :: dit = 1
  logical, parameter :: tderivs = .False.
  integer, parameter :: nord_dt = 4


end module options
