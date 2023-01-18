#include "tests.h"
#include <bhr/autodiff_2nd.h>

#include <cmath>
#include <cassert>

namespace bhr {

#define CHECK(a, b) { \
    if (!(a)) { \
      std::cerr << "\nFailed:\n" #a "\n" << b << '\n'; \
      return false; \
    } \
  }
bool test_autodiff_2nd(void) {
  for (int k = 0; k < 100; ++k) {
    double a = random_double(-100.0, 100.0);
    double b = random_double(-1.0, 1.0);
    double c = random_double(-1.0, 1.0);
    {
      first_partial_derivatives<double, 2> X(a, b, c);
      first_partial_derivatives<double, 2> Y = X * X;
      first_partial_derivatives<double, 2> Z = sqr(X);

      for (int i = 0; i < Y.M_; ++i)
        CHECK(std::abs(Y.d[i] - Z.d[i]) < 1e-9, X << Y << Z);
    }

    {
      second_partial_derivatives<double, 2> X(a, b, c, 0.134, 0.534, 0.43);
      second_partial_derivatives<double, 2> Y = X * X;
      second_partial_derivatives<double, 2> Z = sqr(X);
      for (int i = 0; i < Y.M_; ++i)
        CHECK(std::abs(Y.d[i] - Z.d[i]) < 1e-9, X << Y << Z);
    }
  }
  {
    double x = 10;
    double y = 70;
    second_partial_derivatives<double, 2> X(x, 1, 0, 0, 0, 0);
    second_partial_derivatives<double, 2> Y(y, 0, 1, 0, 0, 0);
    auto XY = X * Y * Y;
    CHECK(std::abs(XY.value() - x * y * y) < 1e-9, XY);
    CHECK(std::abs(XY.first(0) - y * y) < 1e-9, XY);
    CHECK(std::abs(XY.first(1) - 2 * x * y) < 1e-9, XY);
    CHECK(std::abs(XY.second(0, 0) - 0) < 1e-9, XY);
    CHECK(std::abs(XY.second(0, 1) - 2 * y) < 1e-9, XY.second(0, 1));
    CHECK(std::abs(XY.second(1, 1) - 2 * x) < 1e-9, XY);
  }
  {
    double x = 10;
    second_partial_derivatives<double, 1> X(x, 1, 0);
    second_partial_derivatives<double, 1> Y = inverse(X);
    CHECK(std::abs(Y.value() - 1. / x) < 1e-9, Y);
    CHECK(std::abs(Y.first(0) - -1. / (x * x)) < 1e-9, Y);
    CHECK(std::abs(Y.second(0, 0) - 2. / (x * x * x)) < 1e-9, Y);
  }

  {
    second_partial_derivatives<double, 3> X(
        5, .132, .541, 6., .625, .236, .465, .3456, .3465, .15);
        // 10, 1, 0, 0, 0, 0, 0, 0, 0, 0);
        // 5, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    second_partial_derivatives<double, 3> Y = inverse(X);
    second_partial_derivatives<double, 3> Z = X * Y;
    second_partial_derivatives<double, 3> W = X; W *= Y;
    CHECK(std::abs(Z.value() - 1) < 1e-9, Z);
    for (int i = 1; i < Y.M_; ++i)
      CHECK(std::abs(Z.d[i]) < 1e-9, "i=" << i << "  Z=" << Z);
    CHECK(std::abs(W.value() - 1) < 1e-9, Z);
    for (int i = 1; i < Y.M_; ++i)
      CHECK(std::abs(W.d[i]) < 1e-9, "i=" << i << "  W=" << W);
  }

  for (int k = 0; k < 100; ++k) {
    second_partial_derivatives<double, 3> X(
        // 5, .132, .541, 6., .625, .236, .465, .3456, .3465, .15);
        // 5, 1, 0, 0, 0, 0, 0, 0, 0, 0);
        5, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    for (int i = 0; i < X.M_; ++i)
      X.d[i] = random_double(0.1, 1.0);

    second_partial_derivatives<double, 3> Y = exp(log(X));
    for (int i = 0; i < Y.M_; ++i)
      CHECK(std::abs(Y.d[i] - X.d[i]) < 1e-9, "i=" << i << "  Y=" << Y);

    second_partial_derivatives<double, 3> Z = coth(X);
    second_partial_derivatives<double, 3>
        W = (exp(X) + exp(-X)) / (exp(X) - exp(-X));
    for (int i = 0; i < Z.M_; ++i) {
      CHECK(std::abs(Z.d[i] - W.d[i]) < 1e-9,
            "i=" << i << "\nZ=" << Z << "\nW=" << W);
    }

    second_partial_derivatives<double, 3> Y1 = sqrt(sqr(X));
    for (int i = 0; i < Y1.M_; ++i)
      CHECK(std::abs(Y1.d[i] - X.d[i]) < 1e-9, "i=" << i << "  Y1=" << Y1);
  }

  // /* Test mixed. */
  // {
  //   typedef first_partial_derivatives<double, 2> inner_t;
  //   typedef second_partial_derivatives<inner_t, 2> outer_t;
  //   outer_t X(
  //     inner_t(667.434, 6.7549e+21, -6.7549e+21),
  //     inner_t(6.7549e+21, 2.00235e+40, -1.16705e+41),
  //     inner_t(-6.7549e+21, -1.16705e+41, 2.13387e+41),
  //     inner_t(1.00117e+40, -2.14945e+59, -7.16482e+58),
  //     inner_t(-1.16705e+41, -1.43296e+59, 3.48407e+60),
  //     inner_t(1.06694e+41, 1.74203e+60, -4.79621e+60)
  //   );
  //   // 667.434 [6.7549e+21 -6.7549e+21] [6.7549e+21 [2.00235e+40 -1.16705e+41] -6.7549e+21 [-1.16705e+41 2.13387e+41]] [[1.00117e+40 [-2.14945e+59 -7.16482e+58] -1.16705e+41 [-1.43296e+59 3.48407e+60]] [1.06694e+41 [1.74203e+60 -4.79621e+60]] ]
  //   // 667.434 [6.7549e+21 -6.7549e+21] [6.7549e+21 [2.00235e+40 -1.16705e+41] -6.7549e+21 [-1.16705e+41 2.13387e+41]] [[1.00117e+40 [-2.14945e+59 -7.16482e+58] -1.16705e+41 [-1.43296e+59 3.48407e+60]] [1.06694e+41 [1.74203e+60 -4.79621e+60]] ]

  //   outer_t Y = coth(X);
  //   std::cerr << Y << '\n';
  //   CHECK(!isnan(Y), Y);
  // }
  return true;
}

}  // namespace bhr
