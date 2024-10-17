cmake_minimum_required(VERSION 3.16)

# Project
project(n-body-sycl)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

file(MAKE_DIRECTORY images)

set(CMAKE_CXX_COMPILER icpx)

# Check for Intel compiler
if(NOT CMAKE_CXX_COMPILER MATCHES "icpx")
    message(FATAL_ERROR "Intel SYCL requires Intel C++ Compiler (icpx or dpcpp). Please set CMAKE_CXX_COMPILER appropriately.")
endif()

# Find SYCL
find_package(oneDPL REQUIRED)

# SDL3 Setup
# 1. Try system SDL3 package first
find_package(SDL3 QUIET)
if(SDL3_FOUND)
    message(STATUS "Using SDL3 via find_package")
endif()

# 2. Try using a vendored SDL library
if(NOT SDL3_FOUND AND EXISTS "${CMAKE_CURRENT_LIST_DIR}/SDL/CMakeLists.txt")
    add_subdirectory(SDL EXCLUDE_FROM_ALL)
    message(STATUS "Using SDL3 via add_subdirectory")
    set(SDL3_FOUND TRUE)
endif()

# 3. Download SDL, and use that.
if(NOT SDL3_FOUND)
    include(FetchContent)
    set(SDL_SHARED TRUE CACHE BOOL "Build a SDL shared library (if available)")
    set(SDL_STATIC TRUE CACHE BOOL "Build a SDL static library (if available)")
    FetchContent_Declare(
        SDL
        GIT_REPOSITORY https://github.com/libsdl-org/SDL.git
        GIT_TAG main  # Replace this with a particular git tag or git hash
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
    )
    message(STATUS "Using SDL3 via FetchContent")
    FetchContent_MakeAvailable(SDL)
    set_property(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/_deps/sdl-src" PROPERTY EXCLUDE_FROM_ALL TRUE)
endif()

# SDL_image Setup
# 1. Try system SDL_image package first
find_package(SDL3_image QUIET)
if(SDL3_IMAGE_FOUND)
    message(STATUS "Using SDL_image via find_package")
endif()

# 2. Try using a vendored SDL_image library
if(NOT SDL3_IMAGE_FOUND AND EXISTS "${CMAKE_CURRENT_LIST_DIR}/SDL_image/CMakeLists.txt")
    add_subdirectory(SDL_image EXCLUDE_FROM_ALL)
    message(STATUS "Using SDL_image via add_subdirectory")
    set(SDL3_IMAGE_FOUND TRUE)
endif()

# 3. Download SDL_image, and use that.
if(NOT SDL3_IMAGE_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        SDL_image
        GIT_REPOSITORY https://github.com/libsdl-org/SDL_image.git
        GIT_TAG main  # Replace this with a particular git tag or git hash
        GIT_SHALLOW TRUE
        GIT_PROGRESS TRUE
    )
    message(STATUS "Using SDL_image via FetchContent")
    FetchContent_MakeAvailable(SDL_image)
    set_property(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/_deps/sdl_image-src" PROPERTY EXCLUDE_FROM_ALL TRUE)
endif()

# Sources
file(GLOB_RECURSE SOURCES src/*.cpp)
file(GLOB_RECURSE HEADERS include/*.hpp include/*.h)

# Main executable
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})

# Link libraries
target_link_libraries(${PROJECT_NAME} PRIVATE SDL3::SDL3 SDL3_image::SDL3_image)

# C++ 17
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)

# Includes
target_include_directories(${PROJECT_NAME} PUBLIC includes)

# # SYCL flags
# target_compile_options(${PROJECT_NAME} PRIVATE -fsycl)
# target_link_options(${PROJECT_NAME} PRIVATE -fsycl)


# SYCL flags for GPU
target_compile_options(${PROJECT_NAME} PRIVATE 
    -fsycl 
     -fsycl-targets=nvptx64-nvidia-cuda
     -Xsycl-target-backend "--cuda-gpu-arch=sm_75"
)
target_link_options(${PROJECT_NAME} PRIVATE 
    -fsycl 
     -fsycl-targets=nvptx64-nvidia-cuda
     -Xsycl-target-backend "--cuda-gpu-arch=sm_75"
)

