module options
  !
  ! A module to set options for run-time in Mescaline, including:
  !     - Aspects of the simulation, e.g. intitial time, timestep...
  !     - Aspects of the data itself, e.g. spacing
  !     - Specifics of the domains we want to take averages in
  !     - Which output data to write
  !
  !       NOTE: when changing any options; make sure to clean and re-compile
  !
  use, intrinsic :: iso_c_binding, only:c_double
  implicit none
  integer,        parameter :: dp = 8            ! to make floats double precision
  real(c_double), parameter :: pi = 4._dp * atan(1._dp)
  integer, parameter :: clen = 100               ! default length of character variables for all characters

  ! -----------------------------------------------------------------------------
  !
  ! Specify some things about the NR simulation we are analysing
  !    (most of these can be found in the .par file of the simulation)
  !
  ! -----------------------------------------------------------------------------
  real(c_double), parameter :: ainit = 1._dp                ! initial scale factor
  real(c_double), parameter :: dtfac = 0.5                  ! Courant factor: dt = dtfac * dx
  real(c_double), parameter :: tinitial = 0.0               ! initial time set in .par file
  real(c_double), parameter :: rhoinit  = 6.58243402355512e-09 ! initial background density - this is used if it=0 has not been run (will only affect \delta)
  !
  ! Set the gauge. This is set in get_dtalp in analytic_solns.f90
  !   --> allowed: "harmonic" (dt(alp)= -F alp^N K), "synchronous" (dt(alp)=0), "conformal" (dt(alp)=adot, FLRW)
  !
  character(len=clen), parameter :: alpgauge = "harmonic"
  real(c_double), parameter ::     harmonicN = 2             ! if (harmonic) --> dt(alp)= -harmonicF alp^harmonicN K (as in ML_BSSN)
  real(c_double), parameter ::     harmonicF = 1._dp / 3._dp ! these define the evolution of the lapse dt(alp) set in parameter files

  ! -----------------------------------------------------------------------------
  !
  ! Specify some things about the data we are reading in
  !
  ! -----------------------------------------------------------------------------
  ! Spacing (in iterations) between HDF5 files, e.g. how often we have 3D dumps
  integer, parameter ::     dit = 39
  !
  ! if (restart): will ignore first nord_dt-1 files passed in
  !         else: will run on first nord_dt-1 files without time derivs, i.e. some approximations
  logical, parameter :: restart = .True.

  ! -----------------------------------------------------------------------------
  !
  ! Some things about the derivatives we will take
  !
  ! -----------------------------------------------------------------------------
  ! Do 4th order accurate spatial derivatives? default is second order
  logical, parameter ::  fourth = .True.
  ! Do time derivatives? For \Theta and \sigma^2 a linear (in v) approx will be used if false
  !      some things with time derivatives set to zero if this is false, e.g. fluidR.
  !      this is used mainly for running on the it=0 data
  logical, parameter :: tderivs = .True.
  ! The order of accuracy of the time derivatives. Note we read in nord_dt+1 timesteps + use backward diff.
  integer, parameter :: nord_dt = 4

  ! -----------------------------------------------------------------------------
  !
  ! Set up some things defining the averaging + domains
  !     "radius" here is in code units
  !
  ! -----------------------------------------------------------------------------
  !
  ! For a single radius of averaging
  !   --> NOTE: to set up an all-box average; set radius = 0. AND nspheres = 1
  integer,        parameter ::     nspheres = 1           ! number of spheres to average over
  real(c_double), parameter ::       radius = 0.         ! radius of averaging, if single radius
  integer,        parameter :: sphere_rseed = 4693894     ! random seed to draw origins of spheres. So long as nspheres and this stay constant, same origins will be drawn.
  ! Type of domain to average; allowed: "sphere", "cube"
  !      (if cube:, "radius" above is L/2 of cube side and nspheres=ncubes, origins drawn in same way)
  character(len=clen), parameter :: domain_type = "sphere"
  !
  ! For a loop over several averaging radii
  logical,        parameter :: looprad = .False.  ! loop over several radii of averaging?
  real(c_double), parameter ::  radmin = 50       ! if (looprad) starting radius
  real(c_double), parameter ::  radmax = 500      ! if (looprad) finishing radius
  real(c_double), parameter ::    drad = 50       ! if (looprad) radius step

  ! -----------------------------------------------------------------------------
  !
  ! Set up some things related to writing output data (tripreports.f90)
  !
  ! -----------------------------------------------------------------------------
  logical, parameter ::        write3D = .False.    ! write 3D data? Writes many 3D arrays to one HDF5 file for each timestep
  logical, parameter ::     writetilde = .False.    ! write \tilde{} quantities to 3D data? e.g. \tilde{\Theta} = (alp / Gam) * \Theta
  logical, parameter ::    writeomegas = .False.     ! write \Omega_m, \Omega_R, \Omega_Q to ascii files
  logical, parameter ::        writeaD = .False.     ! write a_D^b to ascii file/s
  logical, parameter ::    writehubble = .False.     ! write Hubble parameter to ascii file/s
  logical, parameter ::     writedelta = .False.     ! write delta = \delta\rho/\bar{\rho} (only for nspheres>1, make sure rhoinit is correct if you haven't run on it=0 file)
  !
  ! Constraint violation:
  !     "raw" -- will output H and |M|=sqrt{M^i M_i}
  !     "rel" -- will output the normalised violation w.r.t the 'energy scale', i.e. H/[H] and M/[M]
  logical, parameter ::               rawHM = .False.   ! output only "raw" values of Ham+Mom violation
  logical, parameter ::         rawandrelHM = .False.   ! output both "raw" and "rel" data to separate files
  logical, parameter :: constraints_l1norm  = .False.    ! output the L1 norm of H and M (raw and/or rel controlled as above)
  !
  ! Some additional output options
  logical, parameter :: writeQdterms   = .False.   ! write the average of each term (except <Theta>) in Qd to a separate file


end module options
