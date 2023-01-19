#ifndef TESTS_H
#define TESTS_H

#include <bhr/base.h>
#include <bhr/coordinate.h>
#include <bhr/utility.h>
#include <iostream>

namespace bhr {

/* Misc */
bool test_float_helpers();

/* Coordinate systems */
bool test_spherical(void);
bool test_boyer_lindquist(void);
bool test_boyer_lindquist2(void);

/* Spacetimes */
bool test_flat_spherical(void);
bool test_kerr(void);

/* Integration */
bool test_integrate__RKF45(void);

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




/* Utility */
template <typename Coord> inline void random_vector(
    Coord &vec,
    const typename Coord::value_type &low,
    const typename Coord::value_type &high) {
  for (size_t i = 0; i < sizeof(vec) / sizeof(typename Coord::value_type); ++i)
    vec[i] = random_double(low, high);
}

template <typename Coord> inline void _random_spherical_vector(
    Coord &vec,
    const typename Coord::value_type &low,
    const typename Coord::value_type &high) {
  vec.t = random_double(low, high);
  vec.r = fabs(random_double(low, high));
  vec.theta = random_double(0, M_PI);
  vec.phi = random_double(0, 2 * M_PI);
}

template <typename T> inline void random_vector(
    SphericalVector4<T> &vec, const T &low, const T &high) {
  _random_spherical_vector(vec, low, high);
}
template <typename T> inline void random_vector(
    BoyerLindquistVector4<T> &vec, const T &low, const T &high) {
  _random_spherical_vector(vec, low, high);
}

template <typename T,
          typename ParamFunc,
          typename Func>
inline bool evaluate_and_compare(
    ParamFunc param_func, Func func, const T *expected, int N, T epsilon) {
  for (int i = 0; i < N; ++i) {
    auto param = param_func(i);
    T received = func(param_func(i));
    T error = received / expected[i] - 1;
    if (std::isnan(expected[i])
        || std::isnan(received)
        || std::abs(error) > epsilon) {
      std::cerr << "Param: " << param << '\n';
      std::cerr << "Received: " << received << '\n';
      std::cerr << "Expected: " << expected[i] << '\n';
      std::cerr << "Relative error: " << error << '\n';
      std::cerr << "Maximum tolerable error: " << epsilon << '\n';
      return false;
    }
  }
  return true;
}

template <typename T>
inline bool _compare_eq_rel(const T &received,
                            const T &expected,
                            const T &epsilon,
                            const T &abs_epsilon,
                            const char *file,
                            int line) {
  T error = received / expected - 1;
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

}  // namespace bhr

#endif
