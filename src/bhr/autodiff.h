#ifndef AUTODIFF_H
#define AUTODIFF_H

#include <cmath>
#include <iostream>
#include <utility>

#include <bhr/utility.h>

template<typename T, int N>
struct first_partial_derivatives {
  typedef T value_type;
  typedef first_partial_derivatives This;

  static constexpr int M_ = N + 1;
  T d[M_];

  first_partial_derivatives() {
    for (int i = 0; i <= N; ++i)
      d[i] = T();
  }
  first_partial_derivatives(const T &_value) {
    d[0] = _value;
    for (int i = 1; i <= N; ++i)
      d[i] = T();
  }
  template <typename First, typename Second, typename... Args>
  explicit first_partial_derivatives(
      const First &first, const Second &second, Args ...args)
      : d{T(first), T(second), T(args)...} {}


  explicit inline operator T(void) const {
    return d[0];
  }
  inline const T &value(void) const { return d[0]; }
  inline T &value(void) { return d[0]; }
  inline const T &first(int i) const { return d[1 + i]; }
  inline T &first(int i) { return d[1 + i]; }

  inline This operator-() const {
    first_partial_derivatives result;
    for (int i = 0; i <= N; ++i)
      result.d[i] = -d[i];
    return result;
  }

  friend inline This operator+(const This &A, const This &B) {
    first_partial_derivatives result;
    for (int i = 0; i <= N; ++i)
      result.d[i] = A.d[i] + B.d[i];
    return result;
  }
  friend inline This operator-(const This &A, const This &B) {
    first_partial_derivatives result;
    for (int i = 0; i <= N; ++i)
      result.d[i] = A.d[i] - B.d[i];
    return result;
  }
  friend inline This operator*(const This &A, const This &B) {
    first_partial_derivatives result;
    result.d[0] = A.d[0] * B.d[0];
    for (int i = 1; i <= N; ++i)
      result.d[i] = A.d[0] * B.d[i] + A.d[i] * B.d[0];
    return result;
  }
  friend inline This operator/(const This &A, const This &B) {
    T inv_b0 = inverse(B.d[0]);
    first_partial_derivatives result;
    result.d[0] = A.d[0] * inv_b0;
    for (int i = 1; i <= N; ++i)
      result.d[i] = (A.d[i] * B.d[0] - A.d[0] * B.d[i]) * sqr(inv_b0);
    return result;
  }

  inline This &operator+=(const This &B) {
    for (int i = 0; i <= N; ++i)
      d[i] += B.d[i];
    return *this;
  }
  inline This &operator-=(const This &B) {
    for (int i = 0; i <= N; ++i)
      d[i] -= B.d[i];
    return *this;
  }
  inline This &operator*=(const This &B) {
    for (int i = 1; i <= N; ++i)
      d[i] = d[0] * B.d[i] + d[i] * B.d[0];
    d[0] *= B.d[0];
    return *this;
  }
  inline This &operator/=(const This &B) {
    return *this = (*this / B);  /* ... */
  }

  // TODO: Other operators.
  friend inline bool operator<(const This &A, const This &B) {
    return A.d[0] < B.d[0];
  }
  friend inline bool operator<(const This &A, const T &B) {
    return A.d[0] < B;
  }
  friend inline bool operator>(const This &A, const This &B) {
    return A.d[0] > B.d[0];
  }
  friend inline bool operator>(const This &A, const T &B) {
    return A.d[0] > B;
  }
  friend inline bool operator==(const This &A, const This &B) {
    return A.d[0] == B.d[0];
  }
  friend inline bool operator==(const This &A, const T &B) {
    return A.d[0] == B;
  }

  friend std::ostream& operator<<(std::ostream& stream, const This &A) {
    stream << A.d[0] << " [";
    for (int i = 1; i <= N; ++i)
      stream << A.d[i] << (i == N ? ']' : ' ');
    return stream;
  }

  friend inline This abs(const This &x) {
    return x > 0 ? x : -x;
  }

  friend inline This sqr(const This &x) {
    This result;
    result.d[0] = sqr(x.d[0]);
    for (int i = 1; i <= N; ++i)
      result.d[i] = 2 * x.d[0] * x.d[i];
    return result;
  }

  friend inline This sqrt(const This &x) {
    using std::sqrt;
    This result;
    const T inv = inverse(x.d[0]);
    result.d[0] = sqrt(x.d[0]);
    const T inv_sqr__half = .5 * inv * result.d[0];
    for (int i = 1; i <= N; ++i)
      result.d[i] = inv_sqr__half * x.d[i];
    return result;
  }

  friend inline This sin(const This &x) {
    using std::sin;
    using std::cos;
    This result;
    result.d[0] = sin(x.d[0]);
    T cos_value = cos(x.d[0]);
    for (int i = 1; i <= N; ++i)
      result.d[i] = x.d[i] * cos_value;
    return result;
  }

  friend inline This cos(const This &x) {
    using std::sin;
    using std::cos;
    This result;
    result.d[0] = cos(x.d[0]);
    T sin_value = sin(x.d[0]);
    for (int i = 1; i <= N; ++i)
      result.d[i] = -x.d[i] * sin_value;
    return result;
  }

  friend inline This exp(const This &x) {
    using std::exp;
    This result;
    const T exp_value = exp(x.d[0]);
    result.d[0] = exp_value;
    for (int i = 1; i <= N; ++i)
      result.d[i] = x.d[i] * exp_value;
    return result;
  }

  friend inline This log(const This &x) {
    using std::log;
    This result;
    result.d[0] = log(x.d[0]);
    const T inv = inverse(x.d[0]);
    for (int i = 1; i <= N; ++i)
      result.d[i] = inv * x.d[i];
    return result;
  }

  friend inline This sinh(const This &x) {
    using std::sinh;
    using std::cosh;
    This result;
    result.d[0] = sinh(x.d[0]);
    T cosh_value = cosh(x.d[0]);
    for (int i = 1; i <= N; ++i)
      result.d[i] = x.d[i] * cosh_value;
    return result;
  }

  friend inline This coth(const This &x) {
    using std::sinh;
    This result;
    const T coth_value = coth(x.d[0]);
    // const T eta = 1 - sqr(coth_value);
    const T eta = -sqr(inverse(sinh(x.d[0])));
    result.d[0] = coth_value;
    for (int i = 1; i <= N; ++i)
      result.d[i] = eta * x.d[i];
    return result;
  }

  friend inline bool isnan(const This &x) {
    for (int i = 0; i <= N; ++i)
      if (isnan(x.d[i]))
        return true;
    return false;
  }
};

#endif
