cmake_minimum_required(VERSION 3.19)

project(BlackHoleRenderer LANGUAGES CXX)

set(BUILD_DIR     "${CMAKE_CURRENT_BINARY_DIR}")
set(ROOT_DIR      "${CMAKE_CURRENT_LIST_DIR}")
set(SRC_DIR       "${ROOT_DIR}/src/bhr")

# Use Release mode by default.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

# Extension instantiations.
add_executable(
    bhr
    "${SRC_DIR}/3rd/sha1.cpp"
    "${SRC_DIR}/3rd/specrend.cpp"
    "${SRC_DIR}/colors.cpp"
    "${SRC_DIR}/disk.cpp"
    "${SRC_DIR}/eh_gen.cpp"
    "${SRC_DIR}/floatimage.cpp"
    "${SRC_DIR}/kerr.cpp"
    "${SRC_DIR}/main.cpp"
    "${SRC_DIR}/parameters.cpp"
    "${SRC_DIR}/render.cpp"
    "${SRC_DIR}/settings.cpp"
    "${SRC_DIR}/spectrum.cpp"
    "${SRC_DIR}/texture.cpp"
    "${SRC_DIR}/tga.cpp")

target_compile_features(bhr PUBLIC cxx_std_17)
target_compile_options(bhr PRIVATE -Wall -Wextra -march=native)
target_include_directories(bhr PRIVATE "${ROOT_DIR}/src")

# add_subdirectory(tests)
