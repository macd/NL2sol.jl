# Top level CMakeLists.txt file.  Normally, you will not have to mess
# with this but if you are interested in the test suite and/or coverage
# then you can party hard.  The setup here assumes that this
# will be used in an out of source build.  That is, use the following
# (assuming current directory is <root>/deps/src)
#
#    mkdir build
#    cd build
#    cmake ..
#    make         # this makes only the shared library for nl2sol
#                   and installs it at <root>/deps/usr/lib/libnl2sol.so
#
# to run the test suite
#    make test
# 
# this will only generally show:
#------------------------------------------------------
#Running tests...
#Test project /home/ubuntu/Nl2sol.jl/deps/src/build
#    Start 1: nlmain
#1/1 Test #1: nlmain ...........................   Passed    0.05 sec
#
#100% tests passed, 0 tests failed out of 1
#
#Total Test time (real) =   0.05 sec
#------------------------------------------------------
#
# The detailed results can be found in 
#
#        <root>/deps/src/build/Testing/Temporary/LastTest.log
#
# For a coverage measurement using the same test suite, first remove all contents
# from build directory. Then issue:
#
#    cmake .. -DCMAKE_BUILD_TYPE=COVERAGE
#    make
#    make test
#
#  The coverage files for the nl2sol shared library (the *.gcda and *.gcno) 
#  can now be found in <root>/deps/src/build/nl2sol/CMakeFiles/nlmain.dir
#
#  and you can see the coverage by running gcov in that directory
#
#     gcov *.o
# > ...
# > Lines executed:87.68% of 1737
#
# Coverage for the test driver is in
#
#  <root>/deps/src/build/tests/CMakeFiles/nlmain.dir
#
# NOTE: don't forget to do a clean rebuild of the library after measuring
# the coverage.  Using coverage enabled libraries in production is not a
# good idea.

cmake_minimum_required (VERSION 2.8)
project (nl2sol)
enable_language (Fortran)
enable_testing()


# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

message ("Build type: " ${CMAKE_BUILD_TYPE})

# default installation
get_filename_component (default_prefix ".." ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory; by default it installs in the Nl2sol.jl directory."
      FORCE)
message ("CMAKE_INSTALL_PREFIX: " ${CMAKE_INSTALL_PREFIX})

# The version number.
set (Nl2sol.jl_VERSION_MAJOR 0)
set (Nl2sol.jl_VERSION_MINOR 1)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
message ("Fortran compiler: " ${Fortran_COMPILER_NAME})

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    # gfortran
    set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -O3")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fbounds-check")
    set (CMAKE_Fortran_FLAGS_COVERAGE   "-O2 --coverage -fPIC")
# Don't yet know the appropriate flags for the Intel compiler or g77
# elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
#     # ifort (untested)
#     set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
#     set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
# elseif (Fortran_COMPILER_NAME MATCHES "g77")
#     # g77
#     set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
#     set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
    message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
    message ("No optimized Fortran compiler flags are known, we just try -O2...")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_INSTALL_PREFIX}/usr/lib)

if (CMAKE_BUILD_TYPE MATCHES "COVERAGE")
    message("WARNING: This is a coverage build. Remember to rebuild a clean lib when done")
endif ()

add_subdirectory(nl2sol)
add_subdirectory(tests)

