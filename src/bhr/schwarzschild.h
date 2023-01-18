#ifndef SCHWARZSCHILD_H
#define SCHWARZSCHILD_H

#include <cmath>
#include <bhr/spacetime.h>

namespace bhr {

class SchwarzschildSpacetime : public BlackHoleBase<SchwarzschildSpacetime> {
 public:
#if PREDEFINED_PARAMS
  constexpr static double M = BLACK_HOLE_M;
#else
  double M;

  SchwarzschildSpacetime(double M_) : M(M_) {}
#endif

  PARAMS_FUNC_STATIC_CONSTEXPR
  real_t black_hole_radius(void) PARAMS_FUNC_CONST {
    return (2 * PHY_G / sqr(PHY_c)) * M;
  }


  // function template partial specialization not allowed, therefore using
  // dummy argument
  template <typename T>
  Null coord_system_parameters(const T &/* dummy */) const {
    return Null();
  }

  template <typename T>
  inline SphericalVector4<T> geodesic_acceleration(
      const SphericalVector4<T> &position,
      const SphericalVector4<T> &direction) const {
#if PREDEFINED_PARAMS
    constexpr auto M = BLACK_HOLE_M;
#endif
    const auto r_rs = position.r - (2 * PHY_G / sqr(PHY_c)) * M;
    const auto one_over_r = 1 / position.r;
    const auto tt = 1 / (position.r * r_rs);

    const auto sin_theta = sin(position.theta);
    const auto cos_theta = cos(position.theta);
    const auto cot_theta = cos_theta / sin_theta;

    SphericalVector4<T> result;
    result.t =
      2 * direction.t * direction.r * M * tt;
    result.r =
      sqr(direction.t) * M * r_rs * cube(one_over_r)
      - sqr(direction.r) * M * tt
      - sqr(direction.theta) * r_rs
      - sqr(direction.phi) * r_rs * sqr(sin_theta);
    result.theta =
      2 * direction.r * direction.theta * one_over_r
      - sqr(direction.phi) * sin_theta * cos_theta;
    result.phi =
      2 * direction.r * direction.phi * one_over_r
      + 2 * direction.theta * direction.phi * cot_theta;
    return result;
  }

  template <typename T>
  Matrix4<T> get_metric_ll(const SphericalVector4<T> &position) const {
    Matrix4<T> result = Matrix4<T>();

    const auto rr = sqr(position.r);
    const auto g_tt = 1 - black_hole_radius() / position.r;

    result[0][0] = -g_tt;
    result[1][1] = 1 / g_tt;
    result[2][2] = rr;
    result[3][3] = rr * sqr(sin(position.theta));

    return result;
  }

  template <typename T>
  Matrix4<T> get_metric_uu(const SphericalVector4<T> &position) const {
    Matrix4<T> result = Matrix4<T>();

    const auto one_over_rr = 1 / sqr(position.r);
    const auto g_tt = 1 - black_hole_radius() / position.r;

    result[0][0] = -1 / g_tt;
    result[1][1] = g_tt;
    result[2][2] = one_over_rr;
    result[3][3] = one_over_rr / sqr(sin(position.theta));

    return result;
  }

  template <typename T>
  inline SphericalVector4<T> get_potential_l(
      const SphericalVector4<T> &position_u) const {
    (void)position_u;
    return {0, 0, 0, 0};
  }
};

// http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?2007ragt.meet....1B&amp;data_type=PDF_HIGH&amp;whole_paper=YES&amp;type=PRINTER&amp;filetype=.pdf
class SchwarzschildDipole {
  const SchwarzschildSpacetime &spacetime;
  const double rs;
 public:
  double mu;
  SchwarzschildDipole(const SchwarzschildSpacetime &_spacetime,
                      double surface_B)
      : spacetime(_spacetime), rs(_spacetime.black_hole_radius()),
        mu(1e29 * UNIT_J / UNIT_T) {
    double surface_max_B = _surface_magnetic_field().second;
    mu *= surface_B / surface_max_B;
    _debug_info();
  }

  template <typename T>
  inline std::pair<T, T> _f(const T &r) const {
    using std::log;
    const T x = r / rs;
    const T y = rs / r;
    const T tmplog = log(1 - y);
    return std::make_pair(
        3 * cube(x) * (tmplog + y * (1 + .5 * y)),
        1.5 / rs * (((6 * x - 3) * x - 1) / (x - 1) + 6 * sqr(x) * tmplog)
    );
  }

  template <typename T>
  inline Matrix4<T> get_F_ll(const SphericalVector4<T> &position_u) const {
    using std::sin;
    using std::log;

    const T x = std::max(position_u.r, (T)NEUTRON_STAR_r) / rs;
    const T tmplog = log(1 - inverse(x));

    const T Frp = -3 * mu / sqr(rs) * sqr(sin(position_u.theta)) * (
        (2 * x - 1) / (x - 1) + 2 * x * tmplog);
    const T Ftp = -1.5 * mu / rs * sin(2 * position_u.theta) * (
        1 + 2 * x * (1 + x * tmplog));

    if (debug >= 3) {
      std::cerr << "x=" << x << "  mu=" << mu << "\n";
      std::cerr << "Frp=" << Frp << "\t"
                << "Ftp=" << Ftp << "\n";
      std::cerr << "Frp=" << Frp / UNIT_T / UNIT_m << " Tm\t"
                << "Ftp=" << Ftp / UNIT_T / sqr(UNIT_m) << " Tm^2\n";
    }

    return Matrix4<T>{{
        {0, 0, 0, 0},
        {0, 0, 0, Frp},
        {0, 0, 0, Ftp},
        {0, -Frp, -Ftp, 0}
    }};
  }

  inline double _calc_F(const SphericalVector4<double> &position_u) const {
    Matrix4<double> F_ll = get_F_ll(position_u);
    Matrix4<double> g_uu = spacetime.get_metric_uu(position_u);
    Matrix4<double> F_lu = F_ll * g_uu;
    Matrix4<double> F_uu = g_uu * F_lu;

    double F = 0;
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        F += F_ll[i][j] * F_uu[i][j];
    F /= 4;
    return F;
  }

  inline double __magnetic_field(double r, double theta) const {
    double phi = random_double(0, 2 * M_PI);
    SphericalVector4<double> position_u{0, r, theta, phi};
    double F = _calc_F(position_u);
    return std::sqrt(2 * F);  // B.
  }

  inline std::pair<double, double> _surface_magnetic_field(void) const {
    constexpr int N = 100;
    double sum = 0;
    double max = 0;
    for (int i = 0; i <= N; ++i) {
      double B = __magnetic_field(NEUTRON_STAR_r, random_double(0, M_PI));
      max = std::max(B, max);
      sum += B;
    }

    return std::make_pair(sum / N, max);
  }

  inline void _debug_info(void) const {
    constexpr int N = 21;
    for (int i = 0; i < N; ++i) {
      double theta = (1e-5 + (1 - 2e-5) * i) * M_PI / (N - 1);
      fprintf(stderr, "   theta=%.2lf deg B=%lg T\n",
          theta * 180 / M_PI,
          __magnetic_field(NEUTRON_STAR_r, theta) / UNIT_T);
    }
    for (int i = 0; i < N; ++i) {
      double r = NEUTRON_STAR_r * (1 + 3. * i / (N / 1));
      fprintf(stderr, "   r=%.2lf deg B=%lg T\n",
          r / NEUTRON_STAR_r, __magnetic_field(r, M_PI / 2) / UNIT_T);
    }

    std::pair<double, double> mag = _surface_magnetic_field();
    fprintf(stderr, "Neutron star r=%lg=%lg km\n",
        NEUTRON_STAR_r, NEUTRON_STAR_r / UNIT_km);
    fprintf(stderr, "Surface magnetic field  avg=%lg T\n", mag.first / UNIT_T);
    fprintf(stderr, "Surface magnetic field  max=%lg T\n", mag.second / UNIT_T);
    fprintf(stderr, "mu=%lg=%lg T m^3\nrs=%lg=%lg km\n",
        mu, mu / (UNIT_T * cube(UNIT_m)), rs, rs / UNIT_km);
  }
};

}  // namespace bhr

#endif
