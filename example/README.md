# Working example

This directory contains a working example of how to use libecpint for its core functionality - calculating integrals and derivatives for, in this example, hydrogen iodide in a CC-pVTZ(-PP) basis.

## Build

To build this, make sure you have built libecpint (and its dependencies, pugixml and Faddeeva). Then using cmake in this directory

```
cmake . -Bbuild -DCMAKE_CXX_FLAGS="-I/include/dirs -std=c++14" -DCMAKE_EXE_LINKER_FLAGS="-L/lib/dirs"
```
where the -I and -L flags are used to point to the libecpint headers and library, respectively, if these are not already in your path.

To then build the library and run the test:
```
cd build
make
./example LIBECPINT_SHARE_DIR
```
where the argument points to the share/libecpint directory.

## With Eigen

If you have Eigen3 installed and want to test the example matrix build, simply add the flag
```
-D_WITH_EIGEN
```
to the CMAKE_CXX_FLAGS in the cmake step. 
