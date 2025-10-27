# By default a shared library will be produced, if a static library is desired, use: 
# "cmake -S . -B build -DVECMATH_BUILD_SHARED=OFF"
# "cmake --build build"

cmake_minimum_required(VERSION 3.15)
project(vecmath C)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Bin)

# user option to control library type
option(VECMATH_BUILD_SHARED "Build vecmath as a shared library" ON)

# create the library
if(VECMATH_BUILD_SHARED)
    add_library(vecmath SHARED 
        vecmath_basic_op.h vecmath_basic_op.c
        vecmath_defines.h
        vecmath_mat_op.h vecmath_mat_op.c
        vecmath_quat_op.h vecmath_quat_op.c
        vecmath_ray_op.h vecmath_ray_op.c
        vecmath_types.h
        vecmath_util.h vecmath_util.c
        vecmath_vec_op.h vecmath_vec_op.c
        vecmath.h
    )
    
    # set compile definitions for shared library build
    target_compile_definitions(vecmath PRIVATE VECMATH_BUILD_SHARED VECMATH_EXPORTS)
    
    # set visibility on linux
    if(UNIX AND NOT ANDROID)
        target_compile_options(vecmath PRIVATE -fvisibility=hidden)
    endif()
else()
    add_library(vecmath STATIC 
        vecmath_basic_op.h vecmath_basic_op.c
        vecmath_defines.h
        vecmath_mat_op.h vecmath_mat_op.c
        vecmath_quat_op.h vecmath_quat_op.c
        vecmath_ray_op.h vecmath_ray_op.c
        vecmath_types.h
        vecmath_util.h vecmath_util.c
        vecmath_vec_op.h vecmath_vec_op.c
        vecmath.h
    )
endif()

# add include directory for public headers
target_include_directories(vecmath PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})