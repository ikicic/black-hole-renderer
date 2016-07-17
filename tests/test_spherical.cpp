#include "tests.h"
#if TESTS_ENABLED

bool test_spherical(void) {
  for (int _test = 0; _test < 100; ++_test) {
    CartesianVector4<double> v, v_old, dv, w_cart;
    random_vector(v, 0.1, 2.0);
    random_vector(dv, -1., 1.);
    SphericalVector4<double> w, dw;
    convert_point_and_diff(Null(), v, dv, Null(), w, dw);

    v_old = v;
    mult_add(&v, 0.0001, dv);
    mult_add(&w, 0.0001, dw);

    convert_point(Null(), w, Null(), w_cart);
    if ((v - w_cart).length() > 0.0000001) {
      std::cerr << "Expected: " << v << " Received: " << w_cart << '\n';
      return false;
    }
  }

  return true;
}
#endif
