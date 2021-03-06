# Architecture

This document describes the high-level architecture of libecpint, to help any new contributors find their way around the codebase.

## High-level overview

The code is roughly divided into four chunks:

- the interface,
- the fixed backend,
- the generated backend,
- the code generator

As a general rule, the API will not change. If functionality is added to the library, it should be exposed through ```api.hpp```, following the style of other functions in that interface. Changing any existing parts of the API will break compatibility, so we will only do this if it is *absolutely* necessary.

The other three parts are discussed below.

## Code Map

### Fixed backend

This is found in ```src/lib``` and can be further divided into the following functionalities:

- Quadrature (```gaussquad.cpp, radial_quad.cpp```)
- Analytic angular integrals (```angular.cpp```)
- Containers (```ecp.cpp, gshell.cpp```)
- Mathematical utilities (```bessel.cpp, mathutil.cpp```)

In brief, this is where you should put changes of any of the above kinds. These are all things that are used by the library as part of core functionality, and that do not change depending on how the generated part of the code is initialised.

Each file generally describes one class (for example, the ```ECP``` object in ```ecp.cpp``` or the ```AngularIntegral``` object in ```angular.cpp```) _OR_ one set of functionalities (e.g. adaptive quadrature for radial integrals, in ```radial_quad.cpp```). Please stick to this convention.

As a general rule, any changes in the fixed backend will only affect other parts of the fixed backend, and will not affect that API or generated code. The exception is ```ecpint.cpp``` which describes an ```ECPIntegral``` object, but connects the fixed and generated code together. As such it is also where much of the integral screening takes place.

### Generated backend

This is found in ```src/generated``` and ```src/generated/radial```. It is where any generated code, or templates required for code generation, are located. Currently there are two "part" files, for the generated file ```qgen.cpp``` that ends up in the main library.

The radial folder is where the unrolling of recurrence relations for the primitive radial integrals happens. Any algorithmic changes with respect to the recurrence relations should be found/placed in here.

### Code generator

The generated code is handled by ```src/generate.cpp``` and the utility functions in ```include/generate.hpp```. The former works by finding all the non-zero integral terms for a particular radial integral class using the ```SumTerm``` objects in the header. These are then sorted and made unique, before the radial code is written into the generated backend described above.   

This part is compiled and run as a separate project before the main library. Therefore if any changes you make depend on the generated code, and in particular the ```qgen``` array in ```ECPIntegral```, you should think very carefully about whether it  will be affected by any changes in the code generation. If so, functionality should be added to the generated backend, in the same way that ```qgen.cpp``` is handled. 
