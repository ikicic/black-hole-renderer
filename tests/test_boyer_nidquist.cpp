#include "tests.h"
#if TESTS_ENABLED
#include "../include/schwarzschild.h"
#include "../include/kerr.h"

bool test_boyer_lindquist2(void) {
  // /* Comparison with spherical for a = 0 */
  // for (int _test = 0; _test < 100; ++_test) {
  //   SphericalVector4<double> v, dv, ddv;
  //   BoyerLindquistVector4<double> w, dw, ddw;

  //   random_vector(v, 1.0, 5.0);
  //   random_vector(dv, -5.0, 5.0);
  //   for (int i = 0; i < 4; ++i) {
  //     w[i] = v[i];
  //     dw[i] = dv[i];
  //   }

  //   ddv = SchwarzschildSpacetime(0.05).geodesic_acceleration(v, dv);
  //   ddw = KerrSpacetime(0.05, 0).geodesic_acceleration(w, dw);

  //   for (int i = 0; i < 4; ++i)
  //     if (fabs(ddv[i] - ddw[i]) > 1e-9) {
  //       std::cerr << "Expected: " << ddv << " Received: " << ddw << '\n';
  //       return false;
  //     }
  // }
  return true;
}

bool test_boyer_lindquist(void) {
#if PREDEFINED_PARAMS
  std::cerr << "Test disabled due to predefined params.";
#else
  for (int _test = 0; _test < 100; ++_test) {
    double a = random_double(0, 0.1);
    CartesianVector4<double> v, v_old, dv, w_cart;
    random_vector(v, -10.0, 10.0);
    random_vector(dv, -1., 1.);
    BoyerLindquistVector4<double> w, dw;
    convert_point_and_diff(Null(), v, dv, a, w, dw);

    v_old = v;
    mult_add(&v, 0.0001, dv);
    mult_add(&w, 0.0001, dw);

    convert_point(a, w, Null(), w_cart);
    if ((v - w_cart).length() > 0.000001) {
      std::cerr << "Expected: " << v << " Received: " << w_cart << '\n';
      return false;
    }
  }

#endif
  return true;
}

#endif
