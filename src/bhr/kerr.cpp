#include <bhr/base.h>

#if CHECK_KERR
#include <bhr/kerr.h>
#include <bhr/schwarzschild.h>
#include <bhr/render.h>

void check_kerr(void) {
  constexpr int s = 1;
  FILE *f = fopen(s == 1 ? "output/kerr_check_prograde.csv"
                         : "output/kerr_check_retrograde.csv", "w");
  assert(f != nullptr);
  assert(NEUTRON_STAR_r < _BLACK_HOLE_r_S);

  constexpr double mdot = PHY_G * BLACK_HOLE_M / sqr(PHY_c);
  constexpr double a = BLACK_HOLE_a;
  const double bsc = -a + 6 * s * mdot * std::cos(std::acos(-s * a / mdot) / 3);

  constexpr double dist = 0.9 * MAX_r;
  OrthographicProjectionCamera camera(
      Vector3{{dist, 0, 0}},
      Vector3{{0, 0, 0}},
      Vector3{{0, 0, 1}},
      1.);
  camera.ver_range = 2 * bsc;

  KerrSpacetime spacetime;
  // SchwarzschildSpacetime spacetime;
  Null field;
  constexpr double dlambda = 0.1;
  constexpr int N = 1000;

  for (int i = 0; i < N; ++i) {
    double bb = i / double(N);
    double b = -bsc / (1 - bb);
    // if (b > dist / 25) break;

    // b = (-20 + 40 * (i / double(N))) * UNIT_km;

    FullGeodesicData<
        BasicGeodesicState<BoyerLindquistVector4<double>>,
//        BasicGeodesicState<SphericalVector4<double>>,
        GeodesicExtraBase<
          GeodesicExtra__Steps,
          GeodesicExtra__DeadReason
        >> result;
    integrate_single_geodesic(
        spacetime,
        field,
        &camera,
        [dlambda](auto, auto){ return dlambda; },
        b / bsc,
        0.,
        &result);
    CartesianVector4<double> pos_cart, dir_cart;
    convert_point_and_diff(
        spacetime.coord_system_parameters(result.basic.position),
        result.basic.position,
        result.basic.direction,
        Null(),
        pos_cart,
        dir_cart);
    double phi_cart = atan2(dir_cart[2], dir_cart[1]);
    double diff = M_PI + phi_cart - result.basic.position.phi;
    if (result.extra.dead_reason != DEAD_BLACK_HOLE) {
      while (diff < -1.9 * M_PI) {
        diff += 2 * M_PI;
        phi_cart += 2 * M_PI;
      }
      while (diff > 2.1 * M_PI) {
        diff -= 2 * M_PI;
        phi_cart -= 2 * M_PI;
      }
    }
    std::cerr << result.extra.steps << " "
              << result.extra.dead_reason << "   "
              << bb << "  "
              << b / UNIT_km << "km  "
              << phi_cart * 180 / M_PI << "deg  "
              << diff * 180 / M_PI << "deg  "
              // << phi_cart * 180 / M_PI << "deg  "
              << pos_cart << '\t'
              << dir_cart << '\n';
              // << result.basic.position << '\t'
              // << result.basic.direction << '\n';
    bool dead = result.extra.dead_reason == DEAD_BLACK_HOLE;
    fprintf(f, "%lg,%lg\n", bb, dead ? -1 : phi_cart);
  }
  std::cerr << "bsc=" << bsc / UNIT_km << "km   "
            << "mdot=" << mdot / UNIT_km << "km\n";

  fclose(f);
}
#endif
