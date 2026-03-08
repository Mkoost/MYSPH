# НИКАКИХ компиляторов и project() здесь не ставим — всё наследуется из корня!

find_package(CUDAToolkit REQUIRED)

add_subdirectory(CUDA)

add_library(SPH INTERFACE)
add_library(SPH::SPH ALIAS SPH)

target_include_directories(SPH INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_compile_features(SPH INTERFACE cxx_std_17)
target_compile_definitions(SPH INTERFACE SPH_HAS_CUDA)

target_link_libraries(SPH INTERFACE SPH_CUDA CUDA::cudart)

message(STATUS "SPH library ready! Use: #include <SPH/SPH.h>")
