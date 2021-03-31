# Mescaline: Extracting useful things from Cactus #

### Purpose ###

* Perform post-simulation analysis on hdf5 datafiles produced by the Cactus code

### Dependencies ###
 
* Several modules from [splash](http://users.monash.edu.au/~dprice/splash/) source code
* HDF5 libraries installed
* C and Fortran compilers

### Compiling Mescaline ###

* Mescaline depends on several files from splash. Hence you must have a copy of the *most recent* splash source code, and the environment variable "SPLASH_DIR" set to the directory containing this, e.g.:
```
export SPLASH_DIR=~/splash
```
* Set your HDF5 location, e.g.:
```
export HDF5ROOT=/usr/local
```
* Type make
```
make
```


### Usage ###

* Run Mescaline on your favourite HDF5 files output from Cactus. Can be run on a single file or a series of files (set appropriate options)
```
bin/mescaline *.h5
```
* NOTE: Ensure you ``make clean'' and re-compile if changing any options in src/options.f90



### Assumptions made ###

* Gauge enforced with shift vector \beta^i = 0 throughout
* Uniform grid dx = dy = dz (nx = ny = nz)
* Periodic boundary conditions implemented for all derivatives


### See the User Guide in doc/ for more details ###
