cmake_minimum_required(VERSION 3.20)
project(nbody LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Detect compiler and set flags
if (MSVC)
  add_compile_options(/O2 /arch:AVX2 /fp:fast)
  # OpenMP (LLVM OpenMP is preferred on Windows for newer features)
  find_package(OpenMP)
  if (OpenMP_CXX_FOUND)
    target_compile_options(nbody PRIVATE ${OpenMP_CXX_FLAGS})
  endif()
else()
  add_compile_options(-O3 -march=native -mavx2 -mfma -ffp-contract=fast -fopenmp)
  add_link_options(-fopenmp)
endif()

add_executable(nbody
  src/main.cpp
  src/force_kernel_avx2.cpp
  src/integrators.cpp
)

# On Windows + MSVC, re-add OpenMP to target (CMake needs target first)
if (MSVC)
  find_package(OpenMP)
  if (OpenMP_CXX_FOUND)
    target_link_libraries(nbody PRIVATE OpenMP::OpenMP_CXX)
  endif()
endif()
