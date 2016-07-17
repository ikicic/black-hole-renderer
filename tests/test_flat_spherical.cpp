#include "../include/base.h"
#include "../include/coordinate.h"
#include "../include/flat.h"

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
