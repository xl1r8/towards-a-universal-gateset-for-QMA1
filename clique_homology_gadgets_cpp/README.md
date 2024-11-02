# Clique homology gadgets in C++

The 4-SAT gadgets from the paper are verified in
[src/test4sat.cpp](src/test4sat.cpp).

## Building

Dependencies: [spectra](https://spectralib.org/), [SuiteSparse](https://people.engr.tamu.edu/davis/suitesparse.html), [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) (already included: [googletest](https://github.com/google/googletest), [cliquer](https://users.aalto.fi/~pat/cliquer.html))

Currently, CMake is configured with a custom build of SuiteSparse, so you need to adapt to your setup.
This is because the SuiteSparse package of my distribution did not seem to have multithreading enabled.
We also provide an alternative [build file](vcpkg_CMakeLists.txt) using [vcpkg](https://vcpkg.io/) to load the dependencies.

The code was tested on Arch Linux and macOS.
Note that 64 GB of RAM are necessary when setting `compute_d1_rank = true` for the 4-qubit gadgets.

## License

All of our code is available under GPLv2 (or later) license.
See [thirdparty](thirdparty) for licenses of dependencies.