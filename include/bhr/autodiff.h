#ifndef AUTODIFF_H
#define AUTODIFF_H

#include <cmath>
#include <iostream>
#include <utility>

#include <bhr/utility.h>

template<typename _T, int _N>
struct first_partial_derivatives {
  typedef _T value_type;
  typedef first_partial_derivatives _This;

  static constexpr int _M = _N + 1;
  _T d[_M];

  first_partial_derivatives() {
    for (int i = 0; i <= _N; ++i)
      d[i] = _T();
  }
  first_partial_derivatives(const _T &_value) {
    d[0] = _value;
    for (int i = 1; i <= _N; ++i)
      d[i] = _T();
  }
  template <typename _First, typename _Second, typename... _Args>
  explicit first_partial_derivatives(
      const _First &first, const _Second &second, _Args ...args)
      : d{_T(first), _T(second), _T(args)...} {}


  explicit inline operator _T(void) const {
    return d[0];
  }
  inline const _T &value(void) const { return d[0]; }
  inline _T &value(void) { return d[0]; }
  inline const _T &first(int i) const { return d[1 + i]; }
  inline _T &first(int i) { return d[1 + i]; }

  inline _This operator-() const {
    first_partial_derivatives result;
    for (int i = 0; i <= _N; ++i)
      result.d[i] = -d[i];
    return result;
  }

  friend inline _This operator+(const _This &A, const _This &B) {
    first_partial_derivatives result;
    for (int i = 0; i <= _N; ++i)
      result.d[i] = A.d[i] + B.d[i];
    return result;
  }
  friend inline _This operator-(const _This &A, const _This &B) {
    first_partial_derivatives result;
    for (int i = 0; i <= _N; ++i)
      result.d[i] = A.d[i] - B.d[i];
    return result;
  }
  friend inline _This operator*(const _This &A, const _This &B) {
    first_partial_derivatives result;
    result.d[0] = A.d[0] * B.d[0];
    for (int i = 1; i <= _N; ++i)
      result.d[i] = A.d[0] * B.d[i] + A.d[i] * B.d[0];
    return result;
  }
  friend inline _This operator/(const _This &A, const _This &B) {
    _T inv_b0 = inverse(B.d[0]);
    first_partial_derivatives result;
    result.d[0] = A.d[0] * inv_b0;
    for (int i = 1; i <= _N; ++i)
      result.d[i] = (A.d[i] * B.d[0] - A.d[0] * B.d[i]) * sqr(inv_b0);
    return result;
  }

  inline _This &operator+=(const _This &B) {
    for (int i = 0; i <= _N; ++i)
      d[i] += B.d[i];
    return *this;
  }
  inline _This &operator-=(const _This &B) {
    for (int i = 0; i <= _N; ++i)
      d[i] -= B.d[i];
    return *this;
  }
  inline _This &operator*=(const _This &B) {
    for (int i = 1; i <= _N; ++i)
      d[i] = d[0] * B.d[i] + d[i] * B.d[0];
    d[0] *= B.d[0];
    return *this;
  }
  inline _This &operator/=(const _This &B) {
    return *this = (*this / B);  /* ... */
  }

  // TODO: Other operators.
  friend inline bool operator<(const _This &A, const _This &B) {
    return A.d[0] < B.d[0];
  }
  friend inline bool operator<(const _This &A, const _T &B) {
    return A.d[0] < B;
  }
  friend inline bool operator>(const _This &A, const _This &B) {
    return A.d[0] > B.d[0];
  }
  friend inline bool operator>(const _This &A, const _T &B) {
    return A.d[0] > B;
  }
  friend inline bool operator==(const _This &A, const _This &B) {
    return A.d[0] == B.d[0];
  }
  friend inline bool operator==(const _This &A, const _T &B) {
    return A.d[0] == B;
  }

  friend std::ostream& operator<<(std::ostream& stream, const _This &A) {
    stream << A.d[0] << " [";
    for (int i = 1; i <= _N; ++i)
      stream << A.d[i] << (i == _N ? ']' : ' ');
    return stream;
  }

  friend inline _This abs(const _This &x) {
    return x > 0 ? x : -x;
  }

  friend inline _This sqr(const _This &x) {
    _This result;
    result.d[0] = sqr(x.d[0]);
    for (int i = 1; i <= _N; ++i)
      result.d[i] = 2 * x.d[0] * x.d[i];
    return result;
  }

  friend inline _This sqrt(const _This &x) {
    using std::sqrt;
    _This result;
    const _T inv = inverse(x.d[0]);
    result.d[0] = sqrt(x.d[0]);
    const _T inv_sqr__half = .5 * inv * result.d[0];
    for (int i = 1; i <= _N; ++i)
      result.d[i] = inv_sqr__half * x.d[i];
    return result;
  }

  friend inline _This sin(const _This &x) {
    using std::sin;
    using std::cos;
    _This result;
    result.d[0] = sin(x.d[0]);
    _T cos_value = cos(x.d[0]);
    for (int i = 1; i <= _N; ++i)
      result.d[i] = x.d[i] * cos_value;
    return result;
  }

  friend inline _This cos(const _This &x) {
    using std::sin;
    using std::cos;
    _This result;
    result.d[0] = cos(x.d[0]);
    _T sin_value = sin(x.d[0]);
    for (int i = 1; i <= _N; ++i)
      result.d[i] = -x.d[i] * sin_value;
    return result;
  }

  friend inline _This exp(const _This &x) {
    using std::exp;
    _This result;
    const _T exp_value = exp(x.d[0]);
    result.d[0] = exp_value;
    for (int i = 1; i <= _N; ++i)
      result.d[i] = x.d[i] * exp_value;
    return result;
  }

  friend inline _This log(const _This &x) {
    using std::log;
    _This result;
    result.d[0] = log(x.d[0]);
    const _T inv = inverse(x.d[0]);
    for (int i = 1; i <= _N; ++i)
      result.d[i] = inv * x.d[i];
    return result;
  }

  friend inline _This sinh(const _This &x) {
    using std::sinh;
    using std::cosh;
    _This result;
    result.d[0] = sinh(x.d[0]);
    _T cosh_value = cosh(x.d[0]);
    for (int i = 1; i <= _N; ++i)
      result.d[i] = x.d[i] * cosh_value;
    return result;
  }

  friend inline _This coth(const _This &x) {
    using std::sinh;
    _This result;
    const _T coth_value = coth(x.d[0]);
    // const _T eta = 1 - sqr(coth_value);
    const _T eta = -sqr(inverse(sinh(x.d[0])));
    result.d[0] = coth_value;
    for (int i = 1; i <= _N; ++i)
      result.d[i] = eta * x.d[i];
    return result;
  }

  friend inline bool isnan(const _This &x) {
    for (int i = 0; i <= _N; ++i)
      if (isnan(x.d[i]))
        return true;
    return false;
  }
};

#endif
