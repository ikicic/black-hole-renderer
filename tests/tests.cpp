#include "tests.h"

#include <cstdio>

namespace bhr {

int debug = 1;  // Global variable in the bhr_lib.

#define TEST(func, name) \
  printf("Testing %s... ", name); \
  if (func()) { \
    printf("OK\n"); \
  } else { \
    printf("FAILED\n"); \
    return false; \
  }

static bool test_all(void) {
  TEST(test_float_helpers, "Floating-point helper functions");
  TEST(test_euler_heisenberg_special1, "Euler-Heisenberg Special Function 1");
  TEST(test_euler_heisenberg_special2, "Euler-Heisenberg Special Function 2");
  TEST(test_autodiff_2nd, "Automatic Differentiation 2nd order");
  TEST(test_geodesic_acceleration__magnetic_field,
       "Geodesic acceleration with magnetic fields.");

  TEST(test_spherical, "Spherical");
  TEST(test_boyer_lindquist, "Boyer Lindquist");
  TEST(test_boyer_lindquist2, "Boyer Lindquist 2");
  TEST(test_kerr, "Kerr");
  // TEST(test_flat_spherical, "Flat in spherical");
  TEST(test_integrate__RKF45, "Runge-Kutta-Fehlberg 4(5) Integration");

  // TEST(test_magnetic_only_lagrangian__base, "Magnetic-only Lagrangian base");
  TEST(test_qed_metric_correction_lambda__low_limit, "Lowest order Lagrangian");

  TEST(test_euler_heisenberg, "Euler-Heisenberg Lagrangian");
  // TEST(test_schwarzschild_dipole, "Schwarzschild Dipole");  UNIMPLEMENTED

  return true;
}

}  // namespace bhr

int main() {
  return bhr::test_all() ? 0 : 1;
}
