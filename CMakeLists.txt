cmake_minimum_required(VERSION 3.19)

project(BlackHoleRenderer LANGUAGES CXX)

set(BUILD_DIR     "${CMAKE_CURRENT_BINARY_DIR}")
set(ROOT_DIR      "${CMAKE_CURRENT_LIST_DIR}")
set(SRC_DIR       "${ROOT_DIR}/src/bhr")
set(TEST_DIR      "${ROOT_DIR}/tests")

option(BHR_ENABLE_TESTS "Should tests be compiled?" ON)

# Use Release mode by default.
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

set(cxx_flags -Wall -Wextra -march=native)

# Extension instantiations.
add_library(
    bhr_lib
    STATIC
    "${SRC_DIR}/3rd/sha1.cpp"
    "${SRC_DIR}/3rd/specrend.cpp"
    "${SRC_DIR}/colors.cpp"
    "${SRC_DIR}/disk.cpp"
    "${SRC_DIR}/eh_gen.cpp"
    "${SRC_DIR}/floatimage.cpp"
    "${SRC_DIR}/kerr.cpp"
    "${SRC_DIR}/parameters.cpp"
    "${SRC_DIR}/render.cpp"
    "${SRC_DIR}/settings.cpp"
    "${SRC_DIR}/spectrum.cpp"
    "${SRC_DIR}/texture.cpp"
    "${SRC_DIR}/tga.cpp")

target_compile_features(bhr_lib PUBLIC cxx_std_17)
target_compile_options(bhr_lib PRIVATE ${cxx_flags})
target_include_directories(bhr_lib PUBLIC "${ROOT_DIR}/src")

add_executable(bhr "${SRC_DIR}/main.cpp")
target_compile_options(bhr PRIVATE ${cxx_flags})
target_link_libraries(bhr PRIVATE bhr_lib)

if (BHR_ENABLE_TESTS)
    add_executable(
        bhr_tests
        "${TEST_DIR}/test_autodiff_2nd.cpp"
        "${TEST_DIR}/test_boyer_nidquist.cpp"
        "${TEST_DIR}/test_euler_heisenberg.cpp"
        "${TEST_DIR}/test_flat.cpp"
        "${TEST_DIR}/test_flat_spherical.cpp"
        "${TEST_DIR}/test_float_helpers.cpp"
        "${TEST_DIR}/test_kerr.cpp"
        "${TEST_DIR}/test_qed_lagrangian.cpp"
        "${TEST_DIR}/test_rgf45.cpp"
        "${TEST_DIR}/test_schwarzschild.cpp"
        "${TEST_DIR}/test_spherical.cpp"
        "${TEST_DIR}/test_tensor.cpp"
        "${TEST_DIR}/tests.cpp"
    )
    target_compile_options(bhr_tests PRIVATE ${cxx_flags} -Wno-unknown-pragmas -Wfloat-conversion -Werror)
    target_include_directories(bhr_tests PRIVATE "${ROOT_DIR}/src")
    target_link_libraries(bhr_tests PRIVATE bhr_lib)

    add_test(NAME bhr_tests COMMAND bhr_tests)
endif()
