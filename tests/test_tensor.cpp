#include "tests.h"
#include <bhr/euler_heisenberg.h>
#include <bhr/field.h>
#include <bhr/flat.h>
#include <bhr/qed_lagrangian.h>
#include <bhr/schwarzschild.h>
#include <bhr/tensor.h>

template <typename _Coord,
          typename _Spacetime,
          typename PotentialLFunc,
          typename FLLFunc>
static bool _test_geoacc(
    const _Spacetime &spacetime,
    PotentialLFunc potential_l_func,
    FLLFunc F_ll_func,
    double epsilon,
    bool compare_empty) {
  // auto F_ll_func = [&](auto position_u) {
  //   typedef typename decltype(position_u)::value_type _T;
  //   typedef first_partial_derivatives<_T, 4> fpds;
  //   CartesianVector4<fpds> ad_position_u{
  //     // fpds(position_u[0], (_T)1, (_T)0, (_T)0, (_T)0),
  //     // fpds(position_u[1], (_T)0, (_T)1, (_T)0, (_T)0),
  //     // fpds(position_u[2], (_T)0, (_T)0, (_T)1, (_T)0),
  //     // fpds(position_u[3], (_T)0, (_T)0, (_T)0, (_T)1),
  //     fpds(position_u[0]),
  //     fpds(position_u[1]),
  //     fpds(position_u[2]),
  //     fpds(position_u[3])
  //   };
  //   ad_position_u[0].first(0) = (_T)1;
  //   ad_position_u[1].first(1) = (_T)1;
  //   ad_position_u[2].first(2) = (_T)1;
  //   ad_position_u[3].first(3) = (_T)1;
  //   const CartesianVector4<fpds> A_l = potential_l_func(ad_position_u);
  //   Matrix4<_T> F_ll;
  //   for (int i = 0; i < 4; ++i)
  //     for (int j = 0; j < 4; ++j)
  //       F_ll[i][j] = A_l[j].first(i) - A_l[i].first(j);
  //   return F_ll;
  // };

  for (int k = 0; k < 100; ++k) {
    _Coord position, direction;
    random_vector(position, 1.0 * NEUTRON_STAR_r, 4.0 * NEUTRON_STAR_r);
    random_vector(direction, -0.1 * NEUTRON_STAR_r, 0.1 * NEUTRON_STAR_r);

    // std::cerr << "\n";
    // std::cerr << "=====================================================\n";
    _Coord result1 = geodesic_acceleration__magnetic_field(
        [&](auto position_u) { return spacetime.get_metric_ll(position_u); },
        [&](auto position_u) { return spacetime.get_metric_uu(position_u); },
        F_ll_func,
        [&](auto F, auto G) { return EH::lagrangian_real(F, G); },
        position,
        direction
    );
    // std::cerr << "-----------------------------------------------------\n";
    _Coord result2 = geodesic_acceleration__magnetic_field__lowest_order(
        [&](auto position_u) { return spacetime.get_metric_ll(position_u); },
        [&](auto position_u) { return spacetime.get_metric_uu(position_u); },
        potential_l_func,
        position,
        direction
    );
    _Coord result3 = !compare_empty
        ? result1
        : spacetime.geodesic_acceleration(position, direction);

    bool ok = true;
    for (int i = 0; i < 4; ++i) {
      ok &= compare_eq_rel(result1[i], result2[i], epsilon, 1e-20);
      ok &= compare_eq_rel(result1[i], result3[i], epsilon, 1e-20);
    }
    if (!ok) {
      std::cerr << "Position  " << position << '\n';
      std::cerr << "Direction  " << direction << '\n';
      std::cerr << "Result 1  " << result1 << '\n';
      std::cerr << "Result 2  " << result2 << '\n';
      std::cerr << "Result 3  " << result3 << '\n';
      return false;
    }
  }

  return true;
}


bool test_geodesic_acceleration__magnetic_field(void) {
  if (NEUTRON_STAR_r == 0.0) {
    printf("Neutron star r=0... SKIPPING! ");
    return true;
  }
  if (!test_flat_magnetic_field()) {
    fprintf(stderr, "Flat Magnetic Field test is required to pass! Aborting!");
    return false;
  }

  /* Compare geodesic accelerations. */
  {
    constexpr double surface_B = 1e3 * UNIT_T;
    FlatDipole field;
    const double factor = surface_B / field._surface_magnetic_field().second;
    auto potential_l_func = [factor, &field](auto position_u) {
      return factor * field.get_potential_l(position_u);
    };
    auto F_ll_func = [factor, &field](auto position_u) {
      return factor * field.get_F_ll(position_u);
    };
    if (!_test_geoacc<CartesianVector4<double>>(
          FlatSpacetime(), potential_l_func, F_ll_func, 1e-12, true)) {
      return false;
    }
  }
  auto zero_potential_l_func = [](auto position_u) {
    return decltype(position_u)();
  };
  auto zero_F_ll_func = [](auto position_u) {
    return Matrix4<typename decltype(position_u)::value_type>();
  };

  if (!_test_geoacc<SphericalVector4<double>>(
        FlatSpacetime(), zero_potential_l_func, zero_F_ll_func, 2e-12, true)) {
    return false;
  }


  /* Schwarzschild metrics */
  SchwarzschildSpacetime spacetime;
  if (!_test_geoacc<SphericalVector4<double>>(
        spacetime, zero_potential_l_func, zero_F_ll_func, 1e-10, true)) {
    return false;
  }

  return true;
}
