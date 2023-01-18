#include "tests.h"
#include <bhr/float_helpers.h>

namespace bhr {

static bool test_constexpr_sqrt() {
  constexpr double expected[] = {
    0.0,
    3.0,
    10.0,
    100.0,
    1.4142135623730950488016887242,
    11.096100666450354441827536898,
  };
  constexpr double computed[] = {
    constexpr_sqrt(0.0),
    constexpr_sqrt(9.0),
    constexpr_sqrt(100.0),
    constexpr_sqrt(10000.0),
    constexpr_sqrt(2.0),
    constexpr_sqrt(123.12345),
  };
  static_assert(sizeof(expected) == sizeof(computed));
  for (size_t i = 0; i < sizeof(expected) / sizeof(expected[0]); ++i) {
    if (!compare_eq_rel(computed[i], expected[i], 1e-15, 1e-16))
      return false;
  }
  return true;
}


bool test_float_helpers() {
  return test_constexpr_sqrt();
}

}  // namespace bhr
