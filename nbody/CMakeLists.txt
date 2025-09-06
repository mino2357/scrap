cmake_minimum_required(VERSION 3.20)
project(nbody LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optimization preset (default|strong|weak)
set(NBODY_OPT_PRESET "default" CACHE STRING "Optimization preset: default|strong|weak")
set_property(CACHE NBODY_OPT_PRESET PROPERTY STRINGS default strong weak)

set(NBODY_CXX_FLAGS "")
set(NBODY_LINK_FLAGS "")

if (MSVC)
  if (NBODY_OPT_PRESET STREQUAL "strong")
    set(NBODY_CXX_FLAGS "/Ox;/arch:AVX2;/fp:fast;/GL")
    set(NBODY_LINK_FLAGS "/LTCG")
  elseif (NBODY_OPT_PRESET STREQUAL "weak")
    set(NBODY_CXX_FLAGS "/Od;/arch:AVX2;/fp:precise")
  else()
    set(NBODY_CXX_FLAGS "/O2;/arch:AVX2;/fp:fast")
  endif()
else()
  if (NBODY_OPT_PRESET STREQUAL "strong")
    set(NBODY_CXX_FLAGS "-Ofast;-march=native;-mavx2;-mfma;-ffast-math;-ffp-contract=fast;-fopenmp;-flto;-fomit-frame-pointer;-funroll-loops")
    set(NBODY_LINK_FLAGS "-fopenmp;-flto")
  elseif (NBODY_OPT_PRESET STREQUAL "weak")
    set(NBODY_CXX_FLAGS "-O0;-g;-march=native;-mavx2;-mfma;-fopenmp;-fno-builtin;-fno-unroll-loops;-fno-omit-frame-pointer")
    set(NBODY_LINK_FLAGS "-fopenmp")
  else()
    set(NBODY_CXX_FLAGS "-O3;-march=native;-mavx2;-mfma;-ffp-contract=fast;-fopenmp")
    set(NBODY_LINK_FLAGS "-fopenmp")
  endif()
endif()

add_executable(nbody
  src/main.cpp
  src/force_kernel_avx2.cpp
  src/integrators.cpp
)

# Apply flags
target_compile_options(nbody PRIVATE ${NBODY_CXX_FLAGS})
target_link_options(nbody PRIVATE ${NBODY_LINK_FLAGS})

message(STATUS "NBODY_OPT_PRESET=${NBODY_OPT_PRESET}")
message(STATUS "NBODY_CXX_FLAGS=${NBODY_CXX_FLAGS}")
message(STATUS "NBODY_LINK_FLAGS=${NBODY_LINK_FLAGS}")

# On Windows + MSVC, re-add OpenMP to target (CMake needs target first)
if (MSVC)
  find_package(OpenMP)
  if (OpenMP_CXX_FOUND)
    target_link_libraries(nbody PRIVATE OpenMP::OpenMP_CXX)
  endif()
endif()
