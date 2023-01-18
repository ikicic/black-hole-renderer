#ifndef COORDINATE_SYSTEMS_H
#define COORDINATE_SYSTEMS_H

#include <bhr/config.h>
#include <bhr/matrix.h>
#include <bhr/mod.h>
#include <bhr/utility.h>

#if PREDEFINED_PARAMS
#include <bhr/parameters.h>
#endif

namespace bhr {

// TODO: Write interface for spherical coordinates?
typedef Vector<real_t, 4> Vector4;
typedef Vector<real_t, 3> Vector3;

template <typename T> struct SphericalPart {
  T theta, phi;

  inline std::pair<double, double> to_txty(void) const {
    // return std::make_pair(spherical.phi, spherical.theta);
    double tx = Mod(phi / (2 * M_PI), 1.);
    if (tx < 0) ++tx;
    double ty = Mod(theta / M_PI, 1.);
    // if (ty < 0) ++ty;
    return std::make_pair(tx, ty);
  }

  inline SphericalPart normalized(void) const {
    double _theta = Mod(theta, 2 * M_PI);
    double _phi;
    if (_theta >= M_PI) {
      _theta -= M_PI;
      _phi = phi + M_PI;
    } else {
      _phi = phi;
    }
    return SphericalPart{_theta, Mod(_phi + M_PI, 2 * M_PI) - M_PI};
  }
};

template <typename T, int N> using CartesianVector = Vector<T, N>;
template <typename T> struct CartesianVector4 : Vector<T, 4> {
  static constexpr bool is_spherical = false;

  typedef T value_type;

  CartesianVector4() = default;
  CartesianVector4(const T &t, const T &x, const T &y, const T &z) {
    (*this)[0] = t;
    (*this)[1] = x;
    (*this)[2] = y;
    (*this)[3] = z;
  }
  CartesianVector4(const Vector<T, 4> &v) : Vector<T, 4>(v) {}

  inline SphericalPart<T> spherical_part(const Null &) const {
    return {
      std::atan2(std::sqrt(sqr((*this)[1]) + sqr((*this)[2])), (*this)[3]),
      std::atan2((*this)[2], (*this)[1])
    };
  }
  inline T get_z(void) const {
    return (*this)[3];
  }
  inline T get_r(void) const {
    return sqrt(sqr((*this)[1]) + sqr((*this)[2]) + sqr((*this)[3]));
  }

  /* Utility specialization BEGIN. */
  friend inline void mult(CartesianVector4 *self, const T &c) {
    for (int i = 0; i < 4; ++i)
      (*self)[i] *= c;
  }

  friend inline void mult_add(
      CartesianVector4 *self, const T &c, const CartesianVector4 &A) {
    for (int i = 0; i < 4; ++i)
      (*self)[i] += c * A[i];
  }

  friend inline void set_and_mult_add(
      CartesianVector4 *self, const CartesianVector4 &A,
      const T &c, const CartesianVector4 &B) {
    for (int i = 0; i < 4; ++i)
      (*self)[i] = A[i] + c * B[i];
  }

  friend inline auto numerical_sqr_distance(
      const CartesianVector4 &A, const CartesianVector4 &B) {
    return sqr(A[0] - B[0]) + sqr(A[1] - B[1]) +
           sqr(A[2] - B[2]) + sqr(A[3] - B[3]);
  }
  /* Utility specialization END. */
};

template <typename T> inline void convert_point(
    const Null &/* cartesian parameters */,
    const CartesianVector4<T> &p,
    const Null &/* cartesian parameters */,
    CartesianVector4<T> &P) {
  P = p;
}

template <typename T> inline void convert_point_and_diff(
    const Null &/* cartesian params */,
    const CartesianVector4<T> &pos_cart,
    const CartesianVector4<T> &diff_cart,
    const Null &/* cartesian params */,
    CartesianVector4<T> &position,
    CartesianVector4<T> &diff) {
  position = pos_cart;
  diff = diff_cart;
}

template <typename T> struct SphericalVector4 {
  static constexpr bool is_spherical = true;
  typedef T value_type;
  T t, r, theta, phi;

  inline SphericalPart<T> spherical_part(const Null &/* dummy */) const {
    return {theta, phi};
  }
  inline T get_z(void) const {
    return r * std::cos(theta);
  }
  inline T get_r(void) const {
    return r;
  }

  inline T &operator[](int k) {
    return ((T *)this)[k];
  }
  inline const T &operator[](int k) const {
    return ((T *)this)[k];
  }

  /* Utility specialization BEGIN. */
  template <typename T2>
  friend inline void mult(SphericalVector4 *self, const T2 &c) {
    self->t *= c;
    self->r *= c;
    self->theta *= c;
    self->phi *= c;
  }

  template <typename T2>
  friend inline void mult_add(SphericalVector4 *self,
      const T2 &c, const SphericalVector4 &A) {
    self->t += c * A.t;
    self->r += c * A.r;
    self->theta += c * A.theta;
    self->phi += c * A.phi;
  }

  template <typename T2>
  friend inline void set_and_mult_add(
      SphericalVector4 *self,
      const SphericalVector4 &A,
      const T2 &c, const SphericalVector4 &B) {
    self->t = A.t + c * B.t;
    self->r = A.r + c * B.r;
    self->theta = A.theta + c * B.theta;
    self->phi = A.phi + c * B.phi;
  }

  friend inline auto numerical_sqr_distance(
      const SphericalVector4 &A, const SphericalVector4 &B) {
    return 0.1 * sqr(A[0] - B[0]) + sqr(A[1] - B[1]) +
           10 * (sqr(A[2] - B[2]) + sqr(A[3] - B[3]));
  }
  /* utility specialization END. */


  friend std::ostream& operator<<(std::ostream& stream, const SphericalVector4 &A) {
    return stream << '(' << A.t << ' ' << A.r << ' '
      << A.theta << ' ' << A.phi << ')';
  }
};

template <typename T> inline void convert_point(
    const Null &/* spherical parameters */,
    const SphericalVector4<T> &p,
    const Null &/* cartesian parameters */,
    CartesianVector4<T> &P) {
  T tmp = p.r * sin(p.theta);
  P[0] = p.t;
  P[1] = tmp * cos(p.phi);
  P[2] = tmp * sin(p.phi);
  P[3] = p.r * cos(p.theta);
}

template <typename T> inline void convert_point(
    const Null &/* cartesian parameters */,
    const CartesianVector4<T> &pos_cart,
    const Null &/* spherical parameters */,
    SphericalVector4<T> &position) {
  auto ss = sqr(pos_cart[1]) + sqr(pos_cart[2]);
  auto zz = sqr(pos_cart[3]);

  position.t = pos_cart[0];
  position.r = sqrt(ss + zz);
  position.theta = atan2(sqrt(ss), pos_cart[3]);
  position.phi = atan2(pos_cart[2], pos_cart[1]);
}

template <typename T> inline void convert_point_and_diff(
    const Null &/* spherical parameters */,
    const SphericalVector4<T> &position,
    const SphericalVector4<T> &diff,
    const Null &/* cartesian parameters */,
    CartesianVector4<T> &pos_cart,
    CartesianVector4<T> &diff_cart) {
  using std::sin;
  using std::cos;
  const auto sint = sin(position.theta);
  const auto cost = cos(position.theta);
  const auto sinp = sin(position.phi);
  const auto cosp = cos(position.phi);
  const auto r = position.r;
  pos_cart[0] = position.t;
  pos_cart[1] = r * sint * cosp;
  pos_cart[2] = r * sint * sinp;
  pos_cart[3] = r * cost;

  Matrix3<T> mat{{
    {sint * cosp, r * cost * cosp, -r * sint * sinp},
    {sint * sinp, r * cost * sinp, r * sint * cosp},
    {cost, -r * sint, 0}}};

  diff_cart[0] = diff[0];
  diff_cart[1] = mat[0][0] * diff[1] + mat[0][1] * diff[2] + mat[0][2] * diff[3];
  diff_cart[2] = mat[1][0] * diff[1] + mat[1][1] * diff[2] + mat[1][2] * diff[3];
  diff_cart[3] = mat[2][0] * diff[1] + mat[2][1] * diff[2] + mat[2][2] * diff[3];
}

template <typename T> inline void convert_point_and_diff(
    const Null &/* cartesian parameters */,
    const CartesianVector4<T> &pos_cart,
    const CartesianVector4<T> &diff_cart,
    const Null &/* spherical parameters */,
    SphericalVector4<T> &position,
    SphericalVector4<T> &diff) {
  using std::sqrt;
  auto ss = sqr(pos_cart[1]) + sqr(pos_cart[2]);
  auto s = sqrt(ss);
  auto zz = sqr(pos_cart[3]);
  auto rr = ss + zz;
  auto r = sqrt(rr);

  position.t = pos_cart[0];
  position.r = r;
  position.theta = atan2(sqrt(ss), pos_cart[3]);
  position.phi = atan2(pos_cart[2], pos_cart[1]);

  auto cost = pos_cart[3] / r;
  auto sint = s / r;
  auto cosp = pos_cart[1] / s;
  auto sinp = pos_cart[2] / s;

  Matrix3<T> mat{{
      {sint * cosp, r * cost * cosp, -r * sint * sinp},
      {sint * sinp, r * cost * sinp, r * sint * cosp},
      {cost, -r * sint, 0}}};
  Matrix3<T> inv = matrix3_inverse(mat);

  diff.t = diff_cart[0];
  diff.r = inv[0][0] * diff_cart[1] + inv[0][1] * diff_cart[2] + inv[0][2] * diff_cart[3];
  diff.theta = inv[1][0] * diff_cart[1] + inv[1][1] * diff_cart[2] + inv[1][2] * diff_cart[3];
  diff.phi = inv[2][0] * diff_cart[1] + inv[2][1] * diff_cart[2] + inv[2][2] * diff_cart[3];
}





struct KerrSpacetimeCoordParams {
#if PREDEFINED_PARAMS
  KerrSpacetimeCoordParams() {}
  KerrSpacetimeCoordParams(const Null &) {}
  constexpr static real_t a = BLACK_HOLE_a;
#else
  real_t a;
#endif
};




template <typename T> struct BoyerLindquistVector4 {
  static constexpr bool is_spherical = true;
  typedef T value_type;
  T t, r, theta, phi;

  inline SphericalPart<T> spherical_part(
      KerrSpacetimeCoordParams /* a */) const {
    return {theta, phi};
  }
  inline T get_z(void) const {
    return r * std::cos(theta);
  }
  inline T get_r(void) const {
    return r;
  }

  inline T &operator[](int k) {
    return ((T *)this)[k];
  }
  inline const T &operator[](int k) const {
    return ((T *)this)[k];
  }

  /* Utility specialization BEGIN. */
  template <typename T2>
  friend inline void mult(BoyerLindquistVector4 *self, const T2 &c) {
    self->t *= c;
    self->r *= c;
    self->theta *= c;
    self->phi *= c;
  }

  template <typename T2>
  friend inline void mult_add(BoyerLindquistVector4 *self,
      const T2 &c, const BoyerLindquistVector4 &A) {
    self->t += c * A.t;
    self->r += c * A.r;
    self->theta += c * A.theta;
    self->phi += c * A.phi;
  }

  template <typename T2>
  friend inline void set_and_mult_add(
      BoyerLindquistVector4 *self,
      const BoyerLindquistVector4 &A,
      const T2 &c, const BoyerLindquistVector4 &B) {
    self->t = A.t + c * B.t;
    self->r = A.r + c * B.r;
    self->theta = A.theta + c * B.theta;
    self->phi = A.phi + c * B.phi;
  }

  friend inline auto numerical_sqr_distance(
      const BoyerLindquistVector4 &A, const BoyerLindquistVector4 &B) {
    return sqr(A[0] - B[0]) + sqr(A[1] - B[1]) +
           sqr(A[2] - B[2]) + sqr(A[3] - B[3]);
  }
  /* utility specialization END. */


  friend std::ostream& operator<<(std::ostream& stream, const BoyerLindquistVector4 &A) {
    return stream << '(' << A.t << ' ' << A.r << ' '
      << A.theta << ' ' << A.phi << ')';
  }
};

template <typename T> inline void convert_point(
    const KerrSpacetimeCoordParams &params,
    const BoyerLindquistVector4<T> &p,
    const Null &/* cartesian parameters */,
    CartesianVector4<T> &P) {
  T tmp = sqrt(sqr(p.r) + sqr(params.a)) * sin(p.theta);
  P[0] = p.t;
  P[1] = tmp * cos(p.phi);
  P[2] = tmp * sin(p.phi);
  P[3] = p.r * cos(p.theta);
}

template <typename T> inline void convert_point_and_diff(
    const KerrSpacetimeCoordParams &params,
    const BoyerLindquistVector4<T> &position,
    const BoyerLindquistVector4<T> &diff,
    const Null &/* cartesian parameters */,
    CartesianVector4<T> &pos_cart,
    CartesianVector4<T> &diff_cart) {
  using std::sin;
  using std::cos;
  using std::sqrt;
  const auto sint = sin(position.theta);
  const auto cost = cos(position.theta);
  const auto sinp = sin(position.phi);
  const auto cosp = cos(position.phi);
  const auto r = position.r;
  const auto rho = sqrt(sqr(r) + sqr(params.a));
  pos_cart[0] = position.t;
  pos_cart[1] = rho * sint * cosp;
  pos_cart[2] = rho * sint * sinp;
  pos_cart[3] = r * cost;

  Matrix3<T> mat{{
    {r / rho * sint * cosp, r * cost * cosp, -r * sint * sinp},
    {r / rho * sint * sinp, r * cost * sinp, r * sint * cosp},
    {cost, -r * sint, 0}}};

  diff_cart[0] = diff[0];
  diff_cart[1] = mat[0][0] * diff[1] + mat[0][1] * diff[2] + mat[0][2] * diff[3];
  diff_cart[2] = mat[1][0] * diff[1] + mat[1][1] * diff[2] + mat[1][2] * diff[3];
  diff_cart[3] = mat[2][0] * diff[1] + mat[2][1] * diff[2] + mat[2][2] * diff[3];
}

template <typename T> inline void convert_point(
    const Null &/* cartesian parameters */,
    const CartesianVector4<T> &pos_cart,
    const KerrSpacetimeCoordParams &params,
    BoyerLindquistVector4<T> &position) {
  using std::sqrt;
  auto ss = sqr(pos_cart[1]) + sqr(pos_cart[2]);
  auto zz = sqr(pos_cart[3]);
  auto aa = sqr(params.a);
  auto RR_aa = ss + zz - aa;
  auto rr = (RR_aa + sqrt(sqr(RR_aa) + 4 * aa * zz)) / 2;
  auto r = sqrt(rr);

  position.t = pos_cart[0];
  position.r = r;
  position.theta = atan2(sqrt(ss) * r, sqrt(rr + aa) * pos_cart[3]);
  position.phi = atan2(pos_cart[2], pos_cart[1]);
}

template <typename T> inline void convert_point_and_diff(
    const Null &/* cartesian parameters */,
    const CartesianVector4<T> &pos_cart,
    const CartesianVector4<T> &diff_cart,
    const KerrSpacetimeCoordParams &params,
    BoyerLindquistVector4<T> &position,
    BoyerLindquistVector4<T> &diff) {
  using std::sqrt;
  auto ss = sqr(pos_cart[1]) + sqr(pos_cart[2]);
  auto zz = sqr(pos_cart[3]);
  auto aa = sqr(params.a);
  auto RR_aa = ss + zz - aa;
  auto rr = (RR_aa + sqrt(sqr(RR_aa) + 4 * aa * zz)) / 2;
  auto r = sqrt(rr);
  auto ra = sqrt(rr + aa);

  position.t = pos_cart[0];
  position.r = r;
  position.theta = atan2(sqrt(ss) * r, sqrt(rr + aa) * pos_cart[3]);
  position.phi = atan2(pos_cart[2], pos_cart[1]);

  auto s = sqrt(ss);
  auto cost = pos_cart[3] / r;
  auto sint = s / ra;
  auto cosp = pos_cart[1] / s;
  auto sinp = pos_cart[2] / s;

  Matrix3<T> mat{{
    {r / ra * sint * cosp, r * cost * cosp, -r * sint * sinp},
    {r / ra * sint * sinp, r * cost * sinp, r * sint * cosp},
    {cost, -r * sint, 0}}};
  Matrix3<T> inv = matrix3_inverse(mat);

  diff.t = diff_cart[0];
  diff.r = inv[0][0] * diff_cart[1] + inv[0][1] * diff_cart[2] + inv[0][2] * diff_cart[3];
  diff.theta = inv[1][0] * diff_cart[1] + inv[1][1] * diff_cart[2] + inv[1][2] * diff_cart[3];
  diff.phi = inv[2][0] * diff_cart[1] + inv[2][1] * diff_cart[2] + inv[2][2] * diff_cart[3];
}

template <typename T>
inline bool isfinite(const CartesianVector4<T> &position) {
  return std::isfinite(position[0])
      && std::isfinite(position[1])
      && std::isfinite(position[2])
      && std::isfinite(position[3]);
}

template <typename T>
inline bool isfinite(const SphericalVector4<T> &position) {
  return std::isfinite(position.t)
      && std::isfinite(position.r)
      && std::isfinite(position.theta)
      && std::isfinite(position.phi);
}

template <typename T>
inline bool isfinite(const BoyerLindquistVector4<T> &position) {
  return std::isfinite(position.t)
      && std::isfinite(position.r)
      && std::isfinite(position.theta)
      && std::isfinite(position.phi);
}

}  // namespace bhr

#endif
