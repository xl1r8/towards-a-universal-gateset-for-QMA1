cmake_minimum_required(VERSION 3.29)
project(clique_homology_gadgets)

set(CMAKE_CXX_STANDARD 20)

find_package(Eigen3 3.4 REQUIRED NO_MODULE)
find_package(Spectra REQUIRED NO_MODULE)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

add_subdirectory(thirdparty/googletest-1.15.2)

# using self-built SuiteSparse, but we can alternatively use a globally installed package
# find_package(SuiteSparse REQUIRED COMPONENTS CHOLMOD SPQR)

# git clone https://github.com/DrTimothyAldenDavis/SuiteSparse
# mkdir build; cd build
# . /opt/intel/oneapi/setvars.sh
# cmake -DBLA_VENDOR=Intel10_64lp ..
# cmake --build . --config Release -j30
set(SUITESPARSE /home/dorian/tmp/SuiteSparse)  # TODO update for your system or use global package

include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR} .
        ${SUITESPARSE}/SPQR/Include
        ${SUITESPARSE}/CHOLMOD/Include
        ${SUITESPARSE}/SuiteSparse_config/
)

set(SRC
        src/simplex.cpp
        src/graph.cpp
        src/gadgets.cpp
        src/homology.cpp
        src/sparsemath.cpp
        thirdparty/cliquer-1.21/cliquer.c
        thirdparty/cliquer-1.21/graph.c
        thirdparty/cliquer-1.21/reorder.c
)

SET(LIBS Eigen3::Eigen Spectra::Spectra
        # SuiteSparse::CHOLMOD SuiteSparse::SPQR
        ${SUITESPARSE}/build/SPQR/libspqr.so ${SUITESPARSE}/build/CHOLMOD/libcholmod.so
)


add_executable(clique_homology_gadgets ${SRC} src/main.cpp)
target_link_libraries(clique_homology_gadgets ${LIBS})


add_executable(clique_homology_gadgets_tests ${SRC} ${TEST_SRC} src/test.cpp src/test4sat.cpp)
target_link_libraries(clique_homology_gadgets_tests ${LIBS} gtest gtest_main)

