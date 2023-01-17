Cmake_minimum_required(VERSION 3.19)

project(BHRTests LANGUAGES CXX)

set(ROOT_DIR  "${CMAKE_CURRENT_LIST_DIR}/..")
set(TEST_DIR  "${ROOT_DIR}/tests")
set(BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}")

# Source files.
add_executable(
    bhr_tests
    "${TEST_DIR}/test_autodiff_2nd.cpp"
    "${TEST_DIR}/test_boyer_nidquist.cpp"
    "${TEST_DIR}/test_euler_heisenberg.cpp"
    "${TEST_DIR}/test_flat.cpp"
    "${TEST_DIR}/test_flat_spherical.cpp"
    "${TEST_DIR}/test_float_helper.cpp"
    "${TEST_DIR}/test_kerr.cpp"
    "${TEST_DIR}/test_qed_lagrangian.cpp"
    "${TEST_DIR}/test_rgf45.cpp"
    "${TEST_DIR}/test_schwarzschild.cpp"
    "${TEST_DIR}/test_spherical.cpp"
    "${TEST_DIR}/test_tensor.cpp"
    "${TEST_DIR}/tests.cpp"
)
target_compile_options(bhr_tests PRIVATE  -Wall -Wextra -Wno-unknown-pragmas -Wfloat-conversion)
target_include_directories(bhr_tests PRIVATE "${ROOT_DIR}/src")

add_test(NAME bhr_tests COMMAND bhr_tests)


