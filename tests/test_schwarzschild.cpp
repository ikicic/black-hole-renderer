#include "tests.h"
#include <bhr/schwarzschild.h>

namespace bhr {

// bool test_schwarzschild_dipole(void) {
//   SchwarzschildSpacetime spacetime;
//   SchwarzschildDipole field(spacetime);
//   for (int k = 0; k < 100; ++k) {
//     SphericalVector4<double> position;
//     double rs = spacetime.black_hole_radius();
//     random_vector(position, rs * 2, rs * 100);
//
//     double F0 = field._calc_F(position);
//     INCORRECT formula (this is for M = 0)
//     double F1 = field.mu * (1 + 3 * sqr(std::cos(position.theta)))
//         / (2 * sqr(cube(position.r)));
//     if (!compare_eq_rel(F0, F1, 1e-9, 1e-20)) {
//       std::cerr << position << '\n';
//       return false;
//     }
//   }
//
//   return true;
// }

}  // namespace bhr
