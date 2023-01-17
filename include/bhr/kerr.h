#ifndef KERR_H
#define KERR_H

#include <bhr/spacetime.h>

void check_kerr(void);

class KerrSpacetime : public BlackHoleBase<KerrSpacetime> {
 public:
#if PREDEFINED_PARAMS
  KerrSpacetime() {
    assert(BLACK_HOLE_a >= 0);
  }
  static constexpr real_t M = BLACK_HOLE_M;
  static constexpr real_t a = BLACK_HOLE_a;
#else
  real_t M;
  real_t a;

  KerrSpacetime(real_t _M, real_t _a) : M(_M), a(_a) {
    assert(a >= 0);
  }
#endif

  PARAMS_FUNC_STATIC_CMATH_CONSTEXPR real_t black_hole_radius(void) PARAMS_FUNC_CONST {
    /* Boyer-Lindquist radius */
    return M < a ? 0 : M + std::sqrt(M * M - a * a);
  }

  // function template partial specialization not allowed, therefore using
  // dummy argument
  template <typename _T> KerrSpacetimeCoordParams coord_system_parameters(
      const BoyerLindquistVector4<_T> & /* dummy */) const {
#if PREDEFINED_PARAMS
    return KerrSpacetimeCoordParams();
#else
    return KerrSpacetimeCoordParams(a);
#endif
  }


  template <typename _T>
  inline Matrix4<_T> keplerian_tetrad_to_coord(
      const BoyerLindquistVector4<_T> &position) const {
    PARAMS_CONSTEXPR auto aa = sqr(a);
    PARAMS_CONSTEXPR auto rs = (2 * PHY_G / sqr(PHY_c)) * M;
    PARAMS_CONSTEXPR auto sqrtM = std::sqrt(M);

    // CONSTANTS MISSING.
    // REFERENCES MISSING.
    const auto r = position.r;
    const auto Omega_K = sqrtM / (r * std::sqrt(r) + a * sqrtM);
    const auto rr = sqr(r);
    const auto rho2 = rr + aa * sqr(std::cos(position.theta));
    const auto Delta = rr - rs * r + aa;
    const auto sin_theta2 = sqr(std::sin(position.theta));
    const auto tmp = rs * r / rho2;
    const auto ttmp = tmp * sin_theta2;

    const auto g_00 = tmp - 1;
    const auto g_03 = -a * ttmp;
    const auto g_33 = (rr + aa + aa * ttmp) * sin_theta2;
    const auto gamma = inv_sqrt(-g_00 - (2 * g_03 + g_33 * Omega_K) * Omega_K);
    const auto lambda = -(g_03 + g_33 * Omega_K) / (g_00 + g_03 * Omega_K);
    const auto tau = inv_sqrt((g_00 * lambda + 2 * g_03) * lambda + g_33);

    return Matrix4<_T>{{
        {gamma, 0, 0, tau * lambda},
        {0, std::sqrt(Delta / rho2), 0, 0},
        {0, 0, -inv_sqrt(rho2), 0},
        {gamma * Omega_K, 0, 0, tau}
    }};
  }


  template <typename _T>
  inline Matrix4<_T> ZAMO_coord_to_tetrad(
      const BoyerLindquistVector4<_T> &position) const {
    PARAMS_CONSTEXPR auto aa = sqr(a);
    PARAMS_CONSTEXPR auto rs = (2 * PHY_G / sqr(PHY_c)) * M;

    const auto r = position.r;
    const auto rr = sqr(r);
    const auto rho2 = rr + aa * sqr(std::cos(position.theta));
    const auto Delta = rr - rs * r + aa;
    const auto sin_theta = std::sin(position.theta);
    const auto Sigma2 = rho2 * (rr + aa) + rs * aa * r * sqr(sin_theta);
    const auto alpha = std::sqrt(rho2 * Delta / Sigma2);
    const auto omega = rs * a * r / Sigma2;

    const auto tmp = std::sqrt(Sigma2 / rho2) * sin_theta;
    return Matrix4<_T>{{
        {alpha, 0, 0, 0},
        {0, std::sqrt(rho2 / Delta), 0, 0},
        {0, 0, std::sqrt(rho2), 0},
        {-tmp * omega, 0, 0, tmp}
    }};
  }

  template <typename _T>
  inline BoyerLindquistVector4<_T> geodesic_acceleration(
      const BoyerLindquistVector4<_T> &position,
      const BoyerLindquistVector4<_T> &direction) const {
    PARAMS_CONSTEXPR auto aa = sqr(a);

    // CONSTANTS MISSING.
    // REFERENCES MISSING.
    const auto r = position.r;
    const auto rr = sqr(r);
    const auto cos_theta = std::cos(position.theta);
    const auto sin_theta = std::sin(position.theta);
    const auto cot_theta = cos_theta / sin_theta;
    const auto sin_cos_theta = sin_theta * cos_theta;
    const auto cos2_theta = sqr(cos_theta);
    const auto sin2_theta = sqr(sin_theta);
    // const auto Sigma = rr + aa * cos2_theta;
    const auto sigma = rr - aa * cos2_theta;
    const auto one_over_Sigma = 1 / (rr + aa * cos2_theta);
    const auto one_over_Sigma2 = sqr(one_over_Sigma);
    const auto Delta = rr - 2 * M * r + aa;
    const auto one_over_Delta = 1 / Delta;
    const auto T = 2 * M * aa * r * sin2_theta;

    BoyerLindquistVector4<_T> result;
    result.t = (
        2 * direction.t * direction.r * M * (rr + aa) * sigma * one_over_Delta
        - 2 * direction.t * direction.theta * 2 * M * aa * r * sin_cos_theta
        + 2 * direction.theta * direction.phi * T * a * sin_cos_theta
        + 2 * direction.r * direction.phi * M * a * sin2_theta * (aa * cos2_theta * (aa - rr) - rr * (aa + 3 * rr)) * one_over_Delta
      ) * one_over_Sigma2;

    result.r = (
        sqr(direction.t) * M * Delta * sigma * one_over_Sigma2
        - 2 * direction.t * direction.phi * Delta * M * a * sin2_theta * sigma * one_over_Sigma2
        + sqr(direction.r) * (r * aa * sin2_theta - M * sigma) * one_over_Delta
        - 2 * direction.r * direction.theta * aa * sin_cos_theta
        - sqr(direction.theta) * r * Delta
        + sqr(direction.phi) * Delta * sin2_theta * (M * aa * sin2_theta * sigma * one_over_Sigma2 - r)
      ) * one_over_Sigma;

    result.theta = (
        sin_cos_theta * (
          -sqr(direction.t) * 2 * M * aa * r * one_over_Sigma2
          + 2 * direction.t * direction.phi * 2 * M * a * r * (rr + aa) * one_over_Sigma2
          + sqr(direction.r) * aa * one_over_Delta
          - sqr(direction.theta) * aa
          - sqr(direction.phi) * (rr + aa + T * (one_over_Sigma + (rr + aa) * one_over_Sigma2))
        )
        + 2 * direction.r * direction.theta * r
      ) * one_over_Sigma;

    result.phi =
      2 * direction.theta * direction.phi * cot_theta * (1 + T * one_over_Sigma2)
      + 2 * direction.t * direction.r * M * a * sigma * one_over_Sigma2 * one_over_Delta
      - 2 * direction.t * direction.theta * 2 * M * a * r * cot_theta * one_over_Sigma2
      + 2 * direction.r * direction.phi * (
          r + M * (
            (sqr(aa * sin_cos_theta) - rr * (rr + aa)) * one_over_Sigma2 - rr * one_over_Sigma
            // sqr(aa * sin_cos_theta) - rr * (Sigma + rr + aa)
          ) // * one_over_Sigma2
        ) * one_over_Delta;

    return result;
  }

  template <typename _T>
  inline Matrix4<_T> get_metric_ll(
      const BoyerLindquistVector4<_T> &position) const {
    PARAMS_CONSTEXPR auto aa = sqr(a);

    const auto rr = sqr(position.r);
    const auto rho2 = rr + aa * sqr(cos(position.theta));
    const auto Delta = rr - 2 * M * position.r + aa;
    const auto sin_theta2 = sqr(sin(position.theta));
    const auto tmp = 2 * M * position.r / rho2;
    const auto ttmp = tmp * sin_theta2;

    // Matrix4<_T> metric = Matrix4<_T>();
    // metric[0][0] = tmp - 1;
    // metric[1][1] = rho2 / Delta;
    // metric[2][2] = rho2;
    // metric[3][3] = (rr + aa + aa * ttmp) * sin_theta2;
    // metric[0][3] = metric[3][0] = -a * ttmp;
    // return metric;
    return Matrix4<_T>{{
        {tmp - 1, 0, 0, -a * ttmp},
        {0, rho2 / Delta, 0, 0},
        {0, 0, rho2, 0},
        {-a * ttmp, 0, 0, (rr + aa + aa * ttmp) * sin_theta2}
    }};
  }

  template <typename _T>
  inline Matrix4<_T> get_metric_uu(
      const BoyerLindquistVector4<_T> &position) const {
    PARAMS_CONSTEXPR auto aa = sqr(a);
    PARAMS_CONSTEXPR auto rs = (2 * PHY_G / sqr(PHY_c)) * M;

    const auto rr = sqr(position.r);
    const auto rho2 = rr + aa * sqr(std::cos(position.theta));
    const auto Delta = rr - 2 * M * position.r + aa;
    const auto sin_theta2 = sqr(std::sin(position.theta));
    const auto A = (rr + aa) * rho2 + rs * aa * position.r * sin_theta2;

    // const auto tmp = 2 * M * position.r / rho2;
    // const auto ttmp = tmp * sin_theta2;
    // Matrix4<_T> metric = Matrix4<_T>();
    // metric[0][0] = tmp - 1;
    // metric[1][1] = Delta / rho2;
    // metric[2][2] = inverse(rho2);
    // metric[3][3] = (rr + aa + aa * ttmp) * sin_theta2;
    // metric[0][3] = metric[3][0] = -a * ttmp;
    // return metric;

    const auto g03 = -rs * a * position.r / (PHY_c * rho2 * Delta);
    return Matrix4<_T>{{
      {-A / (sqr(PHY_c) * rho2 * Delta), 0, 0, g03},
      {0, Delta / rho2, 0, 0},
      {0, 0, inverse(rho2), 0},
      {g03, 0, 0, (inverse(sin_theta2) - aa / Delta) / rho2}
      // {g03, 0, 0, (1 - rs * position.r / rho2) / (Delta * sin_theta2)}
    }};
  }

  template <typename _T>
  inline Christoffel<_T> get_christoffel_ull(
      const BoyerLindquistVector4<_T> &position) const {
    PARAMS_CONSTEXPR auto aa = sqr(a);
    PARAMS_CONSTEXPR auto rs = (2 * PHY_G / sqr(PHY_c)) * M;

    // http://arxiv.org/pdf/0904.4184v3.pdf
    // CONSTANTS_MISSING: PHY_c should be included in t.
    const auto r = position.r;
    const auto rr = sqr(r);
    const auto cos_theta = std::cos(position.theta);
    const auto sin_theta = std::sin(position.theta);
    const auto cot_theta = cos_theta / sin_theta;
    const auto sin_cos_theta = sin_theta * cos_theta;
    const auto cos2_theta = sqr(cos_theta);
    const auto sin2_theta = sqr(sin_theta);
    // const auto Sigma = rr + aa * cos2_theta;
    const auto sigma = rr - aa * cos2_theta;
    const auto one_over_Sigma = 1 / (rr + aa * cos2_theta);
    const auto one_over_Sigma2 = sqr(one_over_Sigma);
    const auto one_over_Sigma3 = one_over_Sigma * one_over_Sigma2;
    const auto Delta = rr - rs * r + aa;
    const auto one_over_Delta = 1 / Delta;

    const auto rtt = .5 * sqr(PHY_c) * rs * Delta * sigma * one_over_Sigma3;
    const auto ttr = .5 * rs * (rr + aa) * sigma * one_over_Sigma2 * one_over_Delta;
    const auto tth = -rs * aa * r * sin_cos_theta * one_over_Sigma2;
    const auto rtp = -.5 * PHY_c * rs * Delta * a * sin2_theta * sigma * one_over_Sigma3;
    const auto rrr = (r * aa * sin2_theta - .5 * rs * sigma) * one_over_Sigma * one_over_Delta;
    const auto rrh = - aa * sin_cos_theta * one_over_Sigma;
    const auto rhh = - r * Delta * one_over_Sigma;
    const auto php = cot_theta * (1 + rs * aa * r * sin2_theta * one_over_Sigma2);
    const auto htt = -sqr(PHY_c) * rs * aa * r * sin_cos_theta * one_over_Sigma3;
    const auto ptr = .5 * PHY_c * rs * a * sigma * one_over_Sigma2 * one_over_Delta;
    const auto pth = -PHY_c * rs * a * r * cot_theta * one_over_Sigma2;
    const auto htp = PHY_c * rs * a * r * (rr + aa) * sin_cos_theta * one_over_Sigma3;
    const auto hrr = aa * sin_cos_theta * one_over_Sigma * one_over_Delta;
    const auto hrh = r * one_over_Sigma;
    const auto hhh = -aa * sin_cos_theta * one_over_Sigma;
    const auto thp = rs * a * aa * r * sin2_theta * sin_cos_theta * one_over_Sigma2 / PHY_c;
    const auto trp = .5 * rs * a * sin2_theta * (aa * cos2_theta * (aa - rr) - rr * (aa + 3 * rr)) * one_over_Sigma2 * one_over_Delta / PHY_c;
    const auto prp = (r + .5 * rs * ((sqr(aa * sin_cos_theta) - rr * (rr + aa)) * one_over_Sigma2 - rr * one_over_Sigma)) * one_over_Delta;
    const auto rpp = Delta * sin2_theta * (-r + .5 * rs * aa * sin2_theta * sigma * one_over_Sigma2) * one_over_Sigma;
    const auto hpp = -sin_cos_theta * (sqr(rr + aa) - aa * Delta * sin2_theta + (rr + aa) * rs * aa * r * sin2_theta * one_over_Sigma) * one_over_Sigma2;

    return Christoffel<_T>{{
      {
        {0, ttr, tth, 0},
        {ttr, 0, 0, trp},
        {tth, 0, 0, 0},
        {0, trp, thp, 0},
      }, {
        {rtt, 0, 0, rtp},
        {0, rrr, rrh, 0},
        {0, rrh, rhh, 0},
        {rtp, 0, 0, rpp},
      }, {
        {htt, 0, 0, htp},
        {0, hrr, hrh, 0},
        {0, hrh, hhh, 0},
        {htp, 0, 0, hpp},
      }, {
        {0, ptr, pth, 0},
        {ptr, 0, 0, prp},
        {pth, 0, 0, php},
        {0, prp, php, 0},
      }
    }};
  }
};


#endif
