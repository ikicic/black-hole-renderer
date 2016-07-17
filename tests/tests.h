#ifndef TESTS_H
#define TESTS_H

#define TESTS_ENABLED 0

#if TESTS_ENABLED
#include "../include/base.h"
#include "../include/coordinate.h"
#include "../include/utility.h"


/* Coordinate systems */
bool test_spherical(void);
bool test_boyer_lindquist(void);
bool test_boyer_lindquist2(void);

/* Spacetimes */
bool test_flat_spherical(void);
bool test_kerr(void);

/* Integration */
bool test_integrate__RGF45(void);

/* Automatic differentiation */
bool test_autodiff_2nd(void);

/* Magnetic field */
bool test_flat_magnetic_field(void);

/* Effective Lagrangian */
// bool test_magnetic_only_lagrangian__base(void);
bool test_qed_metric_correction_lambda__low_limit(void);
bool test_euler_heisenberg_special1(void);
bool test_euler_heisenberg_special2(void);
bool test_euler_heisenberg(void);

bool test_geodesic_acceleration__magnetic_field(void);

bool test_schwarzschild_dipole(void);


bool test_all(void);



/* Utility */
template <typename _Coord> inline void random_vector(
    _Coord &vec,
    const typename _Coord::value_type &low,
    const typename _Coord::value_type &high) {
  for (size_t i = 0; i < sizeof(vec) / sizeof(typename _Coord::value_type); ++i)
    vec[i] = random_double(low, high);
}

template <typename _Coord> inline void _random_spherical_vector(
    _Coord &vec,
    const typename _Coord::value_type &low,
    const typename _Coord::value_type &high) {
  vec.t = random_double(low, high);
  vec.r = fabs(random_double(low, high));
  vec.theta = random_double(0, M_PI);
  vec.phi = random_double(0, 2 * M_PI);
}

template <typename _T> inline void random_vector(
    SphericalVector4<_T> &vec, const _T &low, const _T &high) {
  _random_spherical_vector(vec, low, high);
}
template <typename _T> inline void random_vector(
    BoyerLindquistVector4<_T> &vec, const _T &low, const _T &high) {
  _random_spherical_vector(vec, low, high);
}

template <typename _T>
inline bool evaluate_and_compare(
    auto param_func, auto func, const _T *expected, int N, _T epsilon) {
  for (int i = 0; i < N; ++i) {
    auto param = param_func(i);
    _T received = func(param_func(i));
    _T error = received / expected[i] - 1;
    if (std::isnan(expected[i])
        || std::isnan(received)
        || std::abs(error) > epsilon) {
      std::cerr << "Param: " << param << '\n';
      std::cerr << "Received: " << received << '\n';
      std::cerr << "Expected: " << expected[i] << '\n';
      std::cerr << "Relative error: " << error << '\n';
      return false;
    }
  }
  return true;
}

template <typename _T>
inline bool _compare_eq_rel(const _T &received,
                            const _T &expected,
                            const _T &epsilon,
                            const _T &abs_epsilon,
                            const char *file,
                            int line) {
  _T error = received / expected - 1;
  if ((std::isnan(expected)
      || std::isnan(received)
      || std::abs(error) > epsilon)
      && std::abs(received) > abs_epsilon
      && std::abs(expected) > abs_epsilon) {
    std::cerr << file << ":" << line << '\n';
    std::cerr << "Received: " << received << '\n';
    std::cerr << "Expected: " << expected << '\n';
    std::cerr << "Relative error: " << error << '\n';
    std::cerr << "Accepted error: " << epsilon << '\n';
    return false;
  }
  return true;
}

#define compare_eq_rel(a, b, c, d) \
   (_compare_eq_rel((a), (b), (c), (d), __FILE__, __LINE__))

#endif
#endif
