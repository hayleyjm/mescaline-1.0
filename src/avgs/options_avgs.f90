module options
  !
  ! A module to set options for run-time in Mescaline (averaging)
  !
  !  READ ALL COMMENTS CAREFULLY to
  !    ensure all parameters here are set to desired/correct values before compiling
  !
  !       (when changing any of these you must *clean* and re-compile)
  !
  use, intrinsic :: iso_c_binding, only:c_double
  implicit none

  integer,        parameter :: dp = 8            ! precision of floats
  real(c_double), parameter :: pi = 4._dp * atan(1._dp)
  real(c_double), parameter :: eps_t = 1.e-14    ! Epsilon for floating-pt comparisons of time
                                                 ! (/= builtin epsilon as times may come in from outside at lower precision)
  integer, parameter :: clen = 100               ! default length of character variables for all characters
  !
  ! if (restart) and (tderivs): will ignore first nord_dt-1 files passed in
  !         else: will run on first nord_dt-1 files without time derivs, i.e. some approximations
  logical, parameter :: restart = .False.

  ! -----------------------------------------------------------------------------
  !
  ! Specify some things about the NR simulation we are analysing
  !    (most of these can be found in the .par file of the simulation)
  !
  ! -----------------------------------------------------------------------------
  real(c_double), parameter :: tinitial = 0.5966092901089324     ! initial time
  real(c_double), parameter ::    dtfac = 0.1                     ! Courant factor: dt = dtfac * dx
  !
  ! Specify parameters of the initial FLRW background of the simulation ICs
  !     -- HLinit is currently only used IF you don't average on t=tinit slice (to get initial average rho)
  !      ** to be deprecated
  real(c_double), parameter :: HLinit = 10.55349565835734     ! initial value of HL (Hubble param * box length)

  ! -----------------------------------------------------------------------------
  !
  ! Set the gauge. This is set in get_dtalp in analytic_solns.f90
  !   --> allowed: "harmonic" (dt(alp)= -F alp^N K), "synchronous" (dt(alp)=0), "conformal" (dt(alp)=adot, FLRW)
  !
  ! -----------------------------------------------------------------------------
  character(len=clen), parameter :: alpgauge = "harmonic"
  real(c_double), parameter ::     harmonicN = 2             ! if (harmonic) --> dt(alp)= -harmonicF alp^harmonicN K (as in ML_BSSN)
  real(c_double), parameter ::     harmonicF = 1._dp / 3._dp ! these define the evolution of the lapse dt(alp) set in parameter files

  ! -----------------------------------------------------------------------------
  !
  ! Specify some things about the data we are reading in
  !
  ! -----------------------------------------------------------------------------
  ! Spacing (in iterations) between HDF5 files, e.g. how often we have 3D dumps
  integer, parameter :: dit = 1

  ! -----------------------------------------------------------------------------
  !
  ! Some things about the finite-difference derivatives we will take
  !
  ! -----------------------------------------------------------------------------
  ! Do time derivatives? For \Theta and \sigma^2 a linear (in v) approx will be used if false
  !      some things with time derivatives set to zero if this is false, e.g. fluidR (3R instead).
  !            (this is used as False mainly for running on the it=0 data only)
  logical, parameter :: tderivs = .False.
  ! The order of accuracy of the time derivatives. Note we require nord_dt+1 timesteps + use backward diff.
  integer, parameter :: nord_dt = 4

  ! -----------------------------------------------------------------------------
  !
  ! Set up some things defining the averaging + domains
  !        "radius" here is in code units
  !
  ! -----------------------------------------------------------------------------
  !
  ! For a single radius of averaging
  !   --> To set up an all-box average; set radius = 0. AND nspheres = 1
  !
  integer,        parameter ::     nspheres = 10           ! number of spheres to average over
  real(c_double), parameter ::       radius = 0.1         ! radius of averaging, if single radius
  integer,        parameter :: sphere_rseed = 762361     ! random seed to draw origins of spheres. So long as nspheres and this stay constant, same origins will be drawn.
  ! Type of domain to average; allowed: "sphere", "cube"
  !      (if cube:, "radius" above is L/2 of cube side and nspheres=ncubes, origins drawn in same way)
  character(len=clen), parameter :: domain_type = "sphere"
  !
  ! For a loop over several averaging radii
  logical,        parameter :: looprad = .False.  ! loop over several radii of averaging?
  real(c_double), parameter ::  radmin = 50       ! if (looprad) starting radius
  real(c_double), parameter ::  radmax = 200      ! if (looprad) finishing radius
  real(c_double), parameter ::    drad = 50       ! if (looprad) radius step

  ! -----------------------------------------------------------------------------
  !
  ! Set up some things related to writing output data (in tripreports.f90)
  !
  ! -----------------------------------------------------------------------------
  !
  ! Write 3D data? Writes a set of default 3D output to one HDF5 file for each timestep
  logical, parameter ::        write3D = .True.

  !
  ! Add some special things to your 3D output files
  logical, parameter ::    shear_3Dout = .False.   !   -- add to 3D HDF5 files ALL 10x indep. components of \sigma_{mu,nu} (fluid restframe)
  logical, parameter ::      Rdd_3Dout = .False.   !   -- add to 3D HDF5 files ALL 10x indep. components of 4R_{mu,nu}
  logical, parameter ::      Aij_3Dout = .False.   !   -- add to 3D HDF5 files ALL 6x indep. components of A_{ij} (trace-free part of K_ij)
  logical, parameter ::      EUd_3Dout = .False.   !   -- add to 3D HDF5 files ALL 16x components of E^{mu}_{nu} (Electric Weyl)
  logical, parameter ::  shearUd_3Dout = .False.   !   -- if (shear_3Dout .and. shearUd_3Dout), output \sigma^{mu}_{nu} (16x comps instead of 10x sigmadd)
  logical, parameter ::      gdd_3Dout = .False.   ! also output all (10-3)=7x indep. components of g_{mu,nu}? (-3 for shift=0)

  logical, parameter ::     writetilde = .False.    ! write \tilde{} quantities to 3D data? e.g. \tilde{\Theta} = (alp / Gam) * \Theta
  logical, parameter :: fluid_constraints_3Dout = .False.  ! calculate the constraints (currently only Momentum) in the fluid frame and output 3D data

  !
  ! Write scalar data of averaged quantities:
  logical, parameter ::    writeomegas = .True.     ! write \Omega_m, \Omega_R, \Omega_Q to ascii files
  logical, parameter ::        writeaD = .True.     ! write a_D^b to ascii file/s
  logical, parameter ::    writehubble = .True.     ! write Hubble parameter to ascii file/s
  logical, parameter ::     writedelta = .False.     ! write delta = \delta\rho/\bar{\rho} (only for nspheres>1, make sure HLinit is correct if you haven't run on it=0 file)
  logical, parameter :: writeQdterms   = .False.   ! write the average of each term (except <Theta>) in Qd to a separate file

  !
  ! Constraint violation scalar output:
  !     "raw" -- will output H and |M|=sqrt{M^i M_i}
  !     "rel" -- will output the normalised violation w.r.t the 'energy scale', i.e. H/[H] and M/[M]
  logical, parameter ::               rawHM = .False.   ! output only "raw" values of Ham+Mom violation
  logical, parameter ::         rawandrelHM = .True.   ! output both "raw" and "rel" data to separate files
  logical, parameter :: constraints_l1norm  = .True.    ! output the L1 norm of H and M (raw and/or rel controlled as above)


  ! -----------------------------------------------------------------------------
  !
  ! Set up some things that are ONLY used in analytic_solns code // mainly for test cases
  !
  ! -----------------------------------------------------------------------------
  ! ainit is NOT needed for simulation runs
  real(c_double), parameter ::    ainit = 0.000999000999000999   ! initial scale factor
  ! box_size is only used for get_rho_init calculation
  real(c_double), parameter :: box_size = 1.0   ! length of domain in code units


end module options
