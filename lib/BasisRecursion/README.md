# BasisRecursion

C++ implementation of the basic basis recursion algorithm.

## Dependencies

 * Eigen matrix library (http://eigen.tuxfamily.org)
 * `boost::python` (http://www.boost.org/doc/libs/1_48_0/libs/python)

### Debian

 * Eigen matrix library in packet `libeigen3-dev`
 * `boost::python` in packet `libboost-python1.46-dev`
   (At this time the version 1.46 is the most recent one from sid.)

### Brutus
The following brutus modules must be loaded:

```
module load python/2.7.2
module load boost
```

The Eigen matrix library must be in the current directory (currently, brutus has no global eigen installation)
 * edit the Makefile and uncomment the brutus includes
 * copy an eigen3 directory into the current directory, eg:

```
wget http://bitbucket.org/eigen/eigen/get/3.0.5.tar.gz
tar xvzf 3.0.5.tar.gz
mv eigen-eigen-6e7488e20373 eigen3
```

In order to be able to load the compiled module on Brutus, the library path must be updated:

``export LD_LIBRARY_PATH=/cluster/apps/boost/1_47_0_nompi/x86_64/gcc_4.1.2/lib64:$LD_LIBRARY_PATH``


## Compilation

 * edit the makefile to contain your correct library paths (contains a few predefined options)
 * `make HCS` to compile the HyperCubicShape module
 * `make prog` to compile a small test program using the HyperCubicShape in C++
 * `make evaluate` to compile the EvaluateBasis module
 * `make install` copies the EvaluateBasis.so module to the src/WaveBlockNS folder

## Usage

 * use the class `HagedornWavepacketCpp` to utilize the C++ optimized version of the basis recursion
 * currently, only the `HyperCubicShape` is implemented in C++

