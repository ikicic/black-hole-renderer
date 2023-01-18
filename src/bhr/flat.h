#ifndef FLAT_H
#define FLAT_H

#include <bhr/spacetime.h>

namespace bhr {

class FlatSpacetime : public SpacetimeBase<FlatSpacetime> {
  typedef SpacetimeBase<FlatSpacetime> _Parent;
 public:
  FlatSpacetime() : SpacetimeBase<FlatSpacetime>() {}

  static constexpr real_t black_hole_radius(void) {
    return 0;
  }

  template <typename Vector, typename ...Args>
  inline bool is_in_black_hole(const Vector &position,
                               Args... /* args */) const {
    // TODO: This is a hack.
    return position.get_r() < NEUTRON_STAR_r;
  }

  template <typename Vector>
  inline bool is_in_near_flat(const Vector &position,
                              double flat_measure) const {
    // TODO: This is a hack.
    return position.get_r() > 1 / flat_measure;
  }

  // function template partial specialization not allowed, therefore using
  // dummy argument
  template <typename T>
  Null coord_system_parameters(const T &/* dummy */) const {
    return Null();
  }

  template <typename T> Matrix4<T> _get_flat_metric(void) const {
    return Matrix4<T>{{
        {-1, 0, 0, 0},
        {0, 1, 0, 0},
        {0, 0, 1, 0},
        {0, 0, 0, 1}
    }};
  }

  template <typename T> inline Matrix4<T> get_metric_ll(
      const CartesianVector4<T> &/* position */) const {
    return _get_flat_metric<T>();
  }

  template <typename T> inline Matrix4<T> get_metric_uu(
      const CartesianVector4<T> &/* position */) const {
    return _get_flat_metric<T>();
  }

  template <typename T>
  inline Matrix4<T> get_metric_ll(const SphericalVector4<T> &position) const {
    using std::sin;
    Matrix4<T> metric = _get_flat_metric<T>();

    const auto rr = sqr(position.r);
    metric[2][2] = rr;
    metric[3][3] = rr * sqr(sin(position.theta));
    return metric;
  }

  template <typename T>
  inline Matrix4<T> get_metric_uu(const SphericalVector4<T> &position) const {
    using std::sin;
    Matrix4<T> metric = _get_flat_metric<T>();

    const auto one_over_rr = 1 / sqr(position.r);
    metric[2][2] = one_over_rr;
    metric[3][3] = one_over_rr / sqr(sin(position.theta));
    return metric;
  }

  template <typename T>
  inline CartesianVector4<T> geodesic_acceleration(
      const CartesianVector4<T> &/* position */,
      const CartesianVector4<T> &/* direction */) const {

    return CartesianVector4<T>{0, 0, 0, 0};
  }

  template <typename T>
  inline SphericalVector4<T> geodesic_acceleration(
      const SphericalVector4<T> position,
      const SphericalVector4<T> direction) const {

    SphericalVector4<T> result;
    auto one_over_r = 1 / position.r;
    auto sin_theta = sin(position.theta);
    auto cos_theta = cos(position.theta);

    result.t = 0;
    result.r = position.r * (
        - sqr(direction.theta)
        - sqr(direction.phi * sin_theta)
      );
    result.theta =
      - sqr(direction.phi) * sin_theta * cos_theta
      + 2 * direction.r * direction.theta * one_over_r;
    result.phi = 2 * direction.phi * (
        direction.r * one_over_r
        + direction.theta / tan(position.theta)
      );
    return result;
  }

  template <typename Position>
  inline auto get_christoffel_ull(const Position &position) const {
    // The specialization for CartesianVector4 "erases" the general function.
    return (static_cast<const _Parent *>(this))->get_christoffel_ull(position);
  }

  template <typename T>
  inline Christoffel<T> get_christoffel_ull(
      const CartesianVector4<T> &) const {
    return Christoffel<T>();  // Christoffel symbols are trivial.
  }

  template <typename T>
  inline Matrix4<T> keplerian_tetrad_to_coord(
      const SphericalVector4<T> &position) const {
    // PARAMS_CONSTEXPR auto M = BLACK_HOLE_M;
    // PARAMS_CONSTEXPR auto rs = (2 * PHY_G / sqr(PHY_c)) * M;
    constexpr auto M = 0.00001;
    constexpr auto rs = 2 * M;

    // CONSTANTS MISSING.
    // REFERENCES MISSING.
    const auto r = position.r;
    const auto Omega_K = std::sqrt(M / cube(r));
    const auto rr = sqr(r);
    const auto Delta = rr - rs * r;
    const auto sin_theta2 = sqr(std::sin(position.theta));

    const auto g_00 = rs / r - 1;
    const auto g_33 = rr * sin_theta2;
    const auto gamma = inv_sqrt(-g_00 - g_33 * sqr(Omega_K));
    const auto lambda = -g_33 * Omega_K / g_00;
    const auto tau = inv_sqrt(g_00 * sqr(lambda) + g_33);

    return Matrix4<T>{{
        {gamma, 0, 0, tau * lambda},
        {0, std::sqrt(Delta) / r, 0, 0},
        {0, 0, -1 / r, 0},
        {gamma * Omega_K, 0, 0, tau}
    }};
  }

  template <typename T>
  inline Matrix4<T> ZAMO_coord_to_tetrad(
      const SphericalVector4<T> &position) const {
    PARAMS_CONSTEXPR auto M = BLACK_HOLE_M;
    PARAMS_CONSTEXPR auto rs = (2 * PHY_G / sqr(PHY_c)) * M;

    const auto r = position.r;
    const auto rr = sqr(r);
    const auto Delta = rr - rs * r;
    const auto sin_theta = std::sin(position.theta);
    const auto alpha = std::sqrt(Delta) / r;

    return Matrix4<T>{{
        {alpha, 0, 0, 0},
        {0, 1 / alpha, 0, 0},
        {0, 0, r, 0},
        {0, 0, 0, r * sin_theta}
    }};
  }
};

}  // namespace bhr

#endif
