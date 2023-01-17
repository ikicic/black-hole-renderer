#include "tests.h"
#if TESTS_ENABLED
#include <bhr/flat.h>
#include <bhr/field.h>

bool test_flat_magnetic_field(void) {
  static int result = -1;
  if (result != -1)
    return result;
  /* Compare the potential A and the tensor F. */
  typedef CartesianVector4<double> _Coord;
  FlatSpacetime spacetime;
  FlatDipole field;
  for (int test = 0; test < 100; ++test) {
    _Coord position_u;
    random_vector(position_u, 0.1, 1.0);

    typedef first_partial_derivatives<double, 4> fpds;
    CartesianVector4<fpds> ad_position_u{
      fpds(position_u[0], 1, 0, 0, 0),
      fpds(position_u[1], 0, 1, 0, 0),
      fpds(position_u[2], 0, 0, 1, 0),
      fpds(position_u[3], 0, 0, 0, 1),
    };
    const CartesianVector4<fpds> A_l = field.get_potential_l(ad_position_u);
    const Matrix4<double> F_ll = field.get_F_ll(position_u);
    Matrix4<double> F2_ll;

    bool ok = true;
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) {
        F2_ll[i][j] = A_l[j].first(i) - A_l[i].first(j);
        ok &= compare_eq_rel(F_ll[i][j], F2_ll[i][j], 1e-9, 0.);
      }

    if (!ok) {
      std::cerr << F_ll << '\n';
      std::cerr << F2_ll << '\n';
      return (result = false);
    }
  }
  return (result = true);
}

// bool test_flat_spherical(void) {
//   FlatSpacetime spacetime;
//
//   Vector4 position_cart{0, 0, 5, 0};
//   Vector4 direction_cart{1, 1, 2, 3};
//
//   SphericalVector4 position, direction, accel;
//   SphericalVector4::point_and_diff_vector(
//       spacetime.coord_system_parameters(SphericalVector4()),
//       position_cart, direction_cart,
//       position, direction);
//
//   geodesic_acceleration(spacetime, position, direction, accel);
//
//   auto to_cart = [](SphericalVector4<double> w) {
//     return CartesianVector4<double>(
//         w.t,
//         w.r * sin(w.theta) * cos(w.phi),
//         w.r * sin(w.theta) * sin(w.phi),
//         w.r * cos(w.theta));
//   };
//
// }
//
#endif
