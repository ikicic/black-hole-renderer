#ifndef REISSNER_NORDSTROM_H
#define REISSNER_NORDSTROM_H

#include "spacetime.h"

class ReissnerNordstromSpacetime
    : public BlackHoleBase<ReissnerNordstromSpacetime> {
 public:
#if PREDEFINED_PARAMS
  ReissnerNordstromSpacetime() {
    __info();
  }
#else
  real_t M;
  real_t Q;

  ReissnerNordstromSpacetime(real_t _M, real_t _Q) : M(_M), Q(_Q) {
    const real_t r_s = _schwarzschild_radius();
    const real_t r_Q_sqr = _charge_radius_sqr();
    if (sqr(r_s) <= 4 * r_Q_sqr) {
      fprintf(stderr, "Charge too large!\n");
      fprintf(stderr, "r_s=%lf r_Q=%lf\n",
          (double)r_s, (double)std::sqrt(r_Q_sqr));
      assert(sqr(r_s) > 4 * r_Q_sqr);
    }
    __info();
  }
#endif
  void __info(void) const {
    const auto r_s = _schwarzschild_radius();
    const auto E = PHY_Coulomb_k * GET_BH_PARAM_THIS(Q) / sqr(r_s);
    fprintf(stderr, "Schwarzschild radius r_s=%lg km\n", r_s / UNIT_km);
    fprintf(stderr, "Electric field at 2r_s=%lg V/m\n", E / (UNIT_V / UNIT_m));
  }

  PARAMS_FUNC_STATIC_CONSTEXPR
  real_t _schwarzschild_radius(void) PARAMS_FUNC_CONST {
    return (2 * PHY_G / sqr(PHY_c)) * GET_BH_PARAM_THIS(M);
  }

  PARAMS_FUNC_STATIC_CONSTEXPR
  real_t _charge_radius_sqr(void) PARAMS_FUNC_CONST {
    // r_Q
    return (PHY_G * PHY_Coulomb_k / sqr(sqr(PHY_c)))
        * sqr(GET_BH_PARAM_THIS(Q));
  }

  PARAMS_FUNC_STATIC_CMATH_CONSTEXPR
  real_t black_hole_radius(void) PARAMS_FUNC_CONST {
    const auto r_s = _schwarzschild_radius();
    const auto r_Q_sqr = _charge_radius_sqr();
    return (r_s + std::sqrt(sqr(r_s) - 4 * r_Q_sqr)) / 2;
  }

  // Function template partial specialization not allowed, therefore using
  // dummy argument.
  template <typename _T>
  Null coord_system_parameters(const _T &/* dummy */) const {
    return Null();
  }

  template <typename _T>
  Matrix4<_T> get_metric_ll(const SphericalVector4<_T> &position) const {
    Matrix4<_T> result = Matrix4<_T>();

    const auto rr = sqr(position.r);
    const auto g_tt = 1 - _schwarzschild_radius() / position.r
        + _charge_radius_sqr() / rr;

    result[0][0] = -g_tt;
    result[1][1] = 1 / g_tt;
    result[2][2] = rr;
    result[3][3] = rr * sqr(sin(position.theta));

    return result;
  }

  template <typename _T>
  Matrix4<_T> get_metric_uu(const SphericalVector4<_T> &position) const {
    Matrix4<_T> result = Matrix4<_T>();

    const auto one_over_rr = 1 / sqr(position.r);
    const auto g_tt = 1 - _schwarzschild_radius() / position.r
        + _charge_radius_sqr() * one_over_rr;

    result[0][0] = -1 / g_tt;
    result[1][1] = g_tt;
    result[2][2] = one_over_rr;
    result[3][3] = one_over_rr / sqr(sin(position.theta));

    return result;
  }

  template <typename _T>
  inline SphericalVector4<_T> get_potential_l(
      const SphericalVector4<_T> &position_u) const {
    (void)position_u;
    // return {0, 0, 0, 0};
    return {
        ((PHY_Coulomb_k / PHY_c) * GET_BH_PARAM_THIS(Q)) / position_u.r,
        0, 0, 0
    };
  }
};


#endif
