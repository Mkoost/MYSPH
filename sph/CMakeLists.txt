cmake_minimum_required(VERSION 3.24)

# ================== КОМПИЛЯТОРЫ — чистый clang (самый стабильный сейчас) ==================
set(CMAKE_CXX_COMPILER "clang++" CACHE STRING "Clang C++ compiler")
set(CMAKE_CUDA_HOST_COMPILER "clang++" CACHE STRING "")
set(CMAKE_CUDA_COMPILER "/usr/local/cuda/bin/nvcc" CACHE STRING "")

set(CMAKE_CUDA_FLAGS_INIT "-allow-unsupported-compiler --expt-relaxed-constexpr" CACHE STRING "" FORCE)
set(CMAKE_CUDA_ARCHITECTURES "native" CACHE STRING "CUDA architectures" FORCE)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(SPH_Processor LANGUAGES CXX CUDA)

message(STATUS "CXX compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "CUDA host compiler: ${CMAKE_CUDA_HOST_COMPILER}")

# ================== Стандарты и режим сборки ==================
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type: Debug или Release" FORCE)
endif()

message(STATUS "=== Build type: ${CMAKE_BUILD_TYPE} ===")

if(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O3 -march=native -ffast-math -DNDEBUG)
else()
    add_compile_options(-O0 -g3 -Wall)
endif()

# ================== Зависимости ==================
find_package(CUDAToolkit REQUIRED)

find_package(VTK REQUIRED COMPONENTS
    CommonCore
    CommonDataModel
    IOXML
)

# ================== Библиотека SPH ==================
add_subdirectory(SPH)

# ================== Основная программа ==================
add_executable(SPHProcessor src/SPHprocessor.cpp)

target_link_libraries(SPHProcessor PRIVATE
    SPH::SPH
    ${VTK_LIBRARIES}
)

vtk_module_autoinit(
    TARGETS SPHProcessor
    MODULES ${VTK_LIBRARIES}
)

target_compile_features(SPHProcessor PRIVATE cxx_std_17)

message(STATUS "SPHProcessor готов! (clang + CUDA 13.1)")
