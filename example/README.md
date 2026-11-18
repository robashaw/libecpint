# Working example

This directory contains a working example of how to use libecpint for its core functionality: calculating the ECP integrals and their gradient for, in this example, hydrogen iodide in a cc-pVDZ(-PP) basis, with the iodine ECP (ECP28MDF) taken from the built-in library. It also shows how to normalise the basis set coefficients, which libecpint uses exactly as given.

It is built and run as one of libecpint's tests (`ctest -R example`).

## Build

To build it on its own, first build and install libecpint. Then, in this directory,

```bash
cmake . -Bbuild -DCMAKE_PREFIX_PATH=/path/to/libecpint/install
cd build
make
./example /path/to/libecpint/install/share/libecpint
```

where `CMAKE_PREFIX_PATH` is the `CMAKE_INSTALL_PREFIX` libecpint was installed to (not needed if that is a standard location), and the argument of `example` is the `share/libecpint` directory.

## With Eigen

If you have Eigen3 installed and want to try mapping the integrals into an Eigen matrix, add

```bash
-DCMAKE_CXX_FLAGS="-D_WITH_EIGEN -I/path/to/eigen3"
```

to the cmake step.
