#ifndef AUTODIFF_2ND_H
#define AUTODIFF_2ND_H

#include <bhr/autodiff.h>

template <typename _T, int _N>
struct second_partial_derivatives {
  typedef _T value_type;
  typedef second_partial_derivatives _This;

  /* f
   * f_x1, f_x2, ..., f_xN,
   * f_x1x1, f_x1x2, ..., f_x1xN,
   *         f_x2x2, ..., f_x2xN,
   *                      f_xNxN
   */
  static constexpr int _M = 1 + _N + _N * (_N + 1) / 2;
  _T d[_M];

  constexpr static inline int _index(int i, int j) {
    /* Assumes i <= j */
    return 1 + _N + _N * i - i * (i + 1) / 2 + j;
  }
  constexpr static inline int _diag_index(int i) {
    return 1 + _N + _N * i - i * (i - 1) / 2;
  }

  second_partial_derivatives() {
    for (int i = 0; i < _M; ++i)
      d[i] = _T();
  }
  second_partial_derivatives(const _T &_value) {
    d[0] = _value;
    for (int i = 1; i < _M; ++i)
      d[i] = _T();
  }
  template <typename _First, typename _Second, typename... _Args>
  explicit second_partial_derivatives(
      const _First &first, const _Second &second, _Args ...args)
      : d{_T(first), _T(second), _T(args)...} {}

  explicit inline operator _T(void) const {
    return d[0];
  }
  inline const _T &value(void) const { return d[0]; }
  inline _T &value(void) { return d[0]; }
  inline const _T &first(int i) const { return d[1 + i]; }
  inline _T &first(int i) { return d[1 + i]; }
  inline _T second(int i, int j) const {
    return i == j ? 2 * d[_diag_index(i)]
                  : d[i < j ? _index(i, j) : _index(j, i)];
  }

  friend inline _This operator+(const _This &A, const _This &B) {
    _This result;
    for (int i = 0; i < _M; ++i)
      result.d[i] = A.d[i] + B.d[i];
    return result;
  }
  friend inline _This operator-(const _This &A, const _This &B) {
    _This result;
    for (int i = 0; i < _M; ++i)
      result.d[i] = A.d[i] - B.d[i];
    return result;
  }
  friend inline _This operator*(const _This &A, const _This &B) {
    _This result;
    result.d[0] = A.d[0] * B.d[0];
    for (int i = 1; i <= _N; ++i)
      result.d[i] = A.d[0] * B.d[i] + A.d[i] * B.d[0];
    for (int i = 0, k = 1 + _N; i < _N; ++i) {
      result.d[k] = A.d[0] * B.d[k] + A.d[1 + i] * B.d[1 + i] + A.d[k] * B.d[0];
      ++k;
      for (int j = i + 1; j < _N; ++j, ++k) {
        result.d[k] = A.d[0] * B.d[k]
                    + A.d[1 + i] * B.d[1 + j]
                    + A.d[1 + j] * B.d[1 + i]
                    + A.d[k] * B.d[0];
      }
    }
    return result;
  }
  friend inline _This inverse(const _This &A) {
    _This result;
    const _T inv = ::inverse(A.d[0]);
    const _T inv_sqr = sqr(inv);
    result.d[0] = inv;
    for (int i = 1; i <= _N; ++i)
      result.d[i] = -inv_sqr * A.d[i];
    for (int i = 0, k = 1 + _N; i < _N; ++i) {
      result.d[k] = inv_sqr * (inv * sqr(A.d[i + 1]) - A.d[k]);
      ++k;
      for (int j = i + 1; j < _N; ++j, ++k)
        result.d[k] = inv_sqr * (2 * inv * A.d[i + 1] * A.d[j + 1] - A.d[k]);
    }
    return result;
  }
  friend inline _This operator/(const _This &A, const _This &B) {
    return A * inverse(B);
  }

  friend inline _This operator+(const _T &A, const _This &B) {
    _This result = B;
    result.d[0] += A;
    return result;
  }
  friend inline _This operator+(const _This &A, const _T &B) {
    _This result = A;
    result.d[0] += B;
    return result;
  }
  friend inline _This operator-(const _This &A, const _T &B) {
    _This result = A;
    result.d[0] -= B;
    return result;
  }
  friend inline _This operator*(const _T &A, const _This &B) {
    _This result;
    for (int i = 0; i < _M; ++i)
      result.d[i] = A * B.d[i];
    return result;
  }
  friend inline _This operator*(const _This &A, const _T &B) {
    _This result;
    for (int i = 0; i < _M; ++i)
      result.d[i] = A.d[i] * B;
    return result;
  }
  friend inline _This operator/(const _This &A, const _T &B) {
    _This result;
    const _T inv = inverse(B);
    for (int i = 0; i < _M; ++i)
      result.d[i] = A.d[i] * inv;
    return result;
  }

  inline _This &operator+=(const _This &B) {
    for (int i = 0; i < _M; ++i)
      d[i] += B.d[i];
    return *this;
  }
  inline _This &operator-=(const _This &B) {
    for (int i = 0; i < _M; ++i)
      d[i] -= B.d[i];
    return *this;
  }
  inline _This &operator*=(const _This &B) {
    for (int i = 0, k = 1 + _N; i < _N; ++i) {
      d[k] = d[0] * B.d[k] + d[1 + i] * B.d[1 + i] + d[k] * B.d[0];
      ++k;
      for (int j = i + 1; j < _N; ++j, ++k) {
        d[k] = d[0] * B.d[k]
             + d[1 + i] * B.d[1 + j]
             + d[1 + j] * B.d[1 + i]
             + d[k] * B.d[0];
      }
    }
    for (int i = 1; i <= _N; ++i)
      d[i] = d[0] * B.d[i] + d[i] * B.d[0];
    d[0] *= B.d[0];
    return *this;
  }
  inline _This &operator/=(const _This &B) {
    return *this = (*this / B);  /* ... */
  }
  inline _This operator-() const {
    _This result;
    for (int i = 0; i < _M; ++i)
      result.d[i] = -d[i];
    return result;
  }

  inline _This &operator+=(const _T &B) {
    d[0] += B;
    return *this;
  }
  inline _This &operator-=(const _T &B) {
    d[0] -= B;
    return *this;
  }
  inline _This &operator*=(const _T &B) {
    for (int i = 0; i < _M; ++i)
      d[i] *= B;
    return *this;
  }

  /* Too lazy to write other operators now... */
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

  friend std::ostream& operator<<(std::ostream& stream, const _This &A) {
    stream << A.d[0] << " [";
    for (int i = 1; i <= _N; ++i)
      stream << A.d[i] << (i == _N ? ']' : ' ');
    stream << " [";
    for (int i = 0, k = 1 + _N; i < _N; ++i) {
      stream << "[";
      for (int j = i; j < _N; ++j, ++k)
        stream << A.d[k] << (j == _N - 1 ? "] " : " ");
    }
    return stream << "]";
  }

  friend inline _This abs(const _This &x) {
    return x > 0 ? x : -x;
  }

  friend inline _This sqr(const _This &x) {
    _This result;
    result.d[0] = sqr(x.d[0]);
    for (int i = 1; i <= _N; ++i)
      result.d[i] = 2 * x.d[0] * x.d[i];
    for (int i = 0, k = 1 + _N; i < _N; ++i) {
      result.d[k] = 2 * x.d[0] * x.d[k] + sqr(x.d[1 + i]);
      ++k;
      for (int j = i + 1; j < _N; ++j, ++k)
        result.d[k] = 2 * x.d[0] * x.d[k] + 2 * x.d[1 + i] * x.d[1 + j];
    }
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
    for (int i = 0, k = 1 + _N; i < _N; ++i) {
      result.d[k] = inv_sqr__half * (x.d[k] - .25 * inv * sqr(x.d[i + 1]));
      ++k;
      for (int j = i + 1; j < _N; ++j, ++k) {
        result.d[k] = inv_sqr__half * (
            x.d[k] - .5 * inv * x.d[i + 1] * x.d[j + 1]);
      }
    }
    return result;
  }

  friend inline _This sin(const _This &x) {
    using std::sin;
    using std::cos;
    _This result;
    const _T sin_value = sin(x.d[0]);
    const _T cos_value = cos(x.d[0]);
    result.d[0] = sin_value;
    for (int i = 1; i <= _N; ++i)
      result.d[i] = x.d[i] * cos_value;
    for (int i = 0, k = 1 + _N; i < _N; ++i) {
      result.d[k] = (-.5 * sin_value) * sqr(x.d[i + 1]) + cos_value * x.d[k];
      ++k;
      for (int j = i + 1; j < _N; ++j, ++k)
        result.d[k] = -sin_value * x.d[i + 1] * x.d[j + 1] + cos_value * x.d[k];
    }
    return result;
  }

  friend inline _This cos(const _This &x) {
    using std::sin;
    using std::cos;
    _This result;
    const _T sin_value = sin(x.d[0]);
    const _T cos_value = cos(x.d[0]);
    result.d[0] = cos_value;
    for (int i = 1; i <= _N; ++i)
      result.d[i] = x.d[i] * -sin_value;
    for (int i = 0, k = 1 + _N; i < _N; ++i) {
      result.d[k] = (-.5 * cos_value) * sqr(x.d[i + 1]) - sin_value * x.d[k];
      ++k;
      for (int j = i + 1; j < _N; ++j, ++k)
        result.d[k] = -cos_value * x.d[i + 1] * x.d[j + 1] - sin_value * x.d[k];
    }
    return result;
  }

  friend inline _This exp(const _This &x) {
    using std::exp;
    _This result;
    const _T exp_value = exp(x.d[0]);
    result.d[0] = exp_value;
    for (int i = 1; i <= _N; ++i)
      result.d[i] = x.d[i] * exp_value;
    for (int i = 0, k = 1 + _N; i < _N; ++i) {
      result.d[k] = exp_value * (.5 * sqr(x.d[i + 1]) + x.d[k]);
      ++k;
      for (int j = i + 1; j < _N; ++j, ++k)
        result.d[k] = exp_value * (x.d[i + 1] * x.d[j + 1] + x.d[k]);
    }
    return result;
  }

  friend inline _This log(const _This &x) {
    using std::log;
    _This result;
    result.d[0] = log(x.d[0]);
    const _T inv = inverse(x.d[0]);
    for (int i = 1; i <= _N; ++i)
      result.d[i] = inv * x.d[i];
    for (int i = 0, k = 1 + _N; i < _N; ++i) {
      result.d[k] = inv * (x.d[k] - .5 * inv * sqr(x.d[i + 1]));
      ++k;
      for (int j = i + 1; j < _N; ++j, ++k)
        result.d[k] = inv * (x.d[k] - inv * x.d[i + 1] * x.d[j + 1]);
    }
    return result;
  }

  friend inline _This coth(const _This &x) {
    using std::sinh;
    _This result;
    const _T coth_value = coth(x.d[0]);
    const _T eta = 1 - sqr(coth_value);
    // const _T eta2 = -sqr(inverse(sinh(x.d[0])));  // This gets NaNs.
    result.d[0] = coth_value;
    for (int i = 1; i <= _N; ++i)
      result.d[i] = x.d[i] * eta;
    for (int i = 0, k = 1 + _N; i < _N; ++i) {
      result.d[k] = eta * (x.d[k] - coth_value * sqr(x.d[i + 1]));
      ++k;
      for (int j = i + 1; j < _N; ++j, ++k)
        result.d[k] = eta * (x.d[k] - 2 * coth_value * x.d[i + 1] * x.d[j + 1]);
    }
    return result;
  }

  friend inline bool isnan(const _This &x) {
    for (int i = 0; i < _M; ++i)
      if (isnan(x.d[i]))
        return true;
    return false;
  }
};

#endif
