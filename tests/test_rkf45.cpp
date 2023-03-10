#include "tests.h"
#include <bhr/integration.h>

namespace bhr {

bool test_integrate__RKF45(void) {
  // Solve f'(x) = 1 + f^2(x), f(0) = 0 --> f(x) = tan(x)
  auto RHS = [](double x) {
    return 1 + x * x;
  };
  constexpr double x = 1.5;
  constexpr double epsilon = 1e-11;
  double result;
  integrate__RKF45(RHS, 0., x, 0.0001, 0.001, 0.1, epsilon, &result);
  CMATH_CONSTEXPR double expected = std::tan(x);
  std::cout << "Difference: " << std::fabs(result - expected) << "   ";
  if (std::fabs(result - expected) > epsilon * (x / 0.0001)) {
    std::cout << "Expected: " << expected << " Received: " << result << '\n';
    return false;
  }

  return true;
}

}  // namespace bhr
