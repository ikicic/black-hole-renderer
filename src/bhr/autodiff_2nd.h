#ifndef AUTODIFF_2ND_H
#define AUTODIFF_2ND_H

#include <bhr/autodiff.h>

template <typename T, int N>
struct second_partial_derivatives {
  typedef T value_type;
  typedef second_partial_derivatives This;

  /* f
   * f_x1, f_x2, ..., f_xN,
   * f_x1x1, f_x1x2, ..., f_x1xN,
   *         f_x2x2, ..., f_x2xN,
   *                      f_xNxN
   */
  static constexpr int M_ = 1 + N + N * (N + 1) / 2;
  T d[M_];

  constexpr static inline int _index(int i, int j) {
    /* Assumes i <= j */
    return 1 + N + N * i - i * (i + 1) / 2 + j;
  }
  constexpr static inline int _diag_index(int i) {
    return 1 + N + N * i - i * (i - 1) / 2;
  }

  second_partial_derivatives() {
    for (int i = 0; i < M_; ++i)
      d[i] = T();
  }
  second_partial_derivatives(const T &_value) {
    d[0] = _value;
    for (int i = 1; i < M_; ++i)
      d[i] = T();
  }
  template <typename First, typename Second, typename... Args>
  explicit second_partial_derivatives(
      const First &first, const Second &second, Args ...args)
      : d{T(first), T(second), T(args)...} {}

  explicit inline operator T(void) const {
    return d[0];
  }
  inline const T &value(void) const { return d[0]; }
  inline T &value(void) { return d[0]; }
  inline const T &first(int i) const { return d[1 + i]; }
  inline T &first(int i) { return d[1 + i]; }
  inline T second(int i, int j) const {
    return i == j ? 2 * d[_diag_index(i)]
                  : d[i < j ? _index(i, j) : _index(j, i)];
  }

  friend inline This operator+(const This &A, const This &B) {
    This result;
    for (int i = 0; i < M_; ++i)
      result.d[i] = A.d[i] + B.d[i];
    return result;
  }
  friend inline This operator-(const This &A, const This &B) {
    This result;
    for (int i = 0; i < M_; ++i)
      result.d[i] = A.d[i] - B.d[i];
    return result;
  }
  friend inline This operator*(const This &A, const This &B) {
    This result;
    result.d[0] = A.d[0] * B.d[0];
    for (int i = 1; i <= N; ++i)
      result.d[i] = A.d[0] * B.d[i] + A.d[i] * B.d[0];
    for (int i = 0, k = 1 + N; i < N; ++i) {
      result.d[k] = A.d[0] * B.d[k] + A.d[1 + i] * B.d[1 + i] + A.d[k] * B.d[0];
      ++k;
      for (int j = i + 1; j < N; ++j, ++k) {
        result.d[k] = A.d[0] * B.d[k]
                    + A.d[1 + i] * B.d[1 + j]
                    + A.d[1 + j] * B.d[1 + i]
                    + A.d[k] * B.d[0];
      }
    }
    return result;
  }
  friend inline This inverse(const This &A) {
    This result;
    const T inv = ::inverse(A.d[0]);
    const T inv_sqr = sqr(inv);
    result.d[0] = inv;
    for (int i = 1; i <= N; ++i)
      result.d[i] = -inv_sqr * A.d[i];
    for (int i = 0, k = 1 + N; i < N; ++i) {
      result.d[k] = inv_sqr * (inv * sqr(A.d[i + 1]) - A.d[k]);
      ++k;
      for (int j = i + 1; j < N; ++j, ++k)
        result.d[k] = inv_sqr * (2 * inv * A.d[i + 1] * A.d[j + 1] - A.d[k]);
    }
    return result;
  }
  friend inline This operator/(const This &A, const This &B) {
    return A * inverse(B);
  }

  friend inline This operator+(const T &A, const This &B) {
    This result = B;
    result.d[0] += A;
    return result;
  }
  friend inline This operator+(const This &A, const T &B) {
    This result = A;
    result.d[0] += B;
    return result;
  }
  friend inline This operator-(const This &A, const T &B) {
    This result = A;
    result.d[0] -= B;
    return result;
  }
  friend inline This operator*(const T &A, const This &B) {
    This result;
    for (int i = 0; i < M_; ++i)
      result.d[i] = A * B.d[i];
    return result;
  }
  friend inline This operator*(const This &A, const T &B) {
    This result;
    for (int i = 0; i < M_; ++i)
      result.d[i] = A.d[i] * B;
    return result;
  }
  friend inline This operator/(const This &A, const T &B) {
    This result;
    const T inv = inverse(B);
    for (int i = 0; i < M_; ++i)
      result.d[i] = A.d[i] * inv;
    return result;
  }

  inline This &operator+=(const This &B) {
    for (int i = 0; i < M_; ++i)
      d[i] += B.d[i];
    return *this;
  }
  inline This &operator-=(const This &B) {
    for (int i = 0; i < M_; ++i)
      d[i] -= B.d[i];
    return *this;
  }
  inline This &operator*=(const This &B) {
    for (int i = 0, k = 1 + N; i < N; ++i) {
      d[k] = d[0] * B.d[k] + d[1 + i] * B.d[1 + i] + d[k] * B.d[0];
      ++k;
      for (int j = i + 1; j < N; ++j, ++k) {
        d[k] = d[0] * B.d[k]
             + d[1 + i] * B.d[1 + j]
             + d[1 + j] * B.d[1 + i]
             + d[k] * B.d[0];
      }
    }
    for (int i = 1; i <= N; ++i)
      d[i] = d[0] * B.d[i] + d[i] * B.d[0];
    d[0] *= B.d[0];
    return *this;
  }
  inline This &operator/=(const This &B) {
    return *this = (*this / B);  /* ... */
  }
  inline This operator-() const {
    This result;
    for (int i = 0; i < M_; ++i)
      result.d[i] = -d[i];
    return result;
  }

  inline This &operator+=(const T &B) {
    d[0] += B;
    return *this;
  }
  inline This &operator-=(const T &B) {
    d[0] -= B;
    return *this;
  }
  inline This &operator*=(const T &B) {
    for (int i = 0; i < M_; ++i)
      d[i] *= B;
    return *this;
  }

  /* Too lazy to write other operators now... */
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

  friend std::ostream& operator<<(std::ostream& stream, const This &A) {
    stream << A.d[0] << " [";
    for (int i = 1; i <= N; ++i)
      stream << A.d[i] << (i == N ? ']' : ' ');
    stream << " [";
    for (int i = 0, k = 1 + N; i < N; ++i) {
      stream << "[";
      for (int j = i; j < N; ++j, ++k)
        stream << A.d[k] << (j == N - 1 ? "] " : " ");
    }
    return stream << "]";
  }

  friend inline This abs(const This &x) {
    return x > 0 ? x : -x;
  }

  friend inline This sqr(const This &x) {
    This result;
    result.d[0] = sqr(x.d[0]);
    for (int i = 1; i <= N; ++i)
      result.d[i] = 2 * x.d[0] * x.d[i];
    for (int i = 0, k = 1 + N; i < N; ++i) {
      result.d[k] = 2 * x.d[0] * x.d[k] + sqr(x.d[1 + i]);
      ++k;
      for (int j = i + 1; j < N; ++j, ++k)
        result.d[k] = 2 * x.d[0] * x.d[k] + 2 * x.d[1 + i] * x.d[1 + j];
    }
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
    for (int i = 0, k = 1 + N; i < N; ++i) {
      result.d[k] = inv_sqr__half * (x.d[k] - .25 * inv * sqr(x.d[i + 1]));
      ++k;
      for (int j = i + 1; j < N; ++j, ++k) {
        result.d[k] = inv_sqr__half * (
            x.d[k] - .5 * inv * x.d[i + 1] * x.d[j + 1]);
      }
    }
    return result;
  }

  friend inline This sin(const This &x) {
    using std::sin;
    using std::cos;
    This result;
    const T sin_value = sin(x.d[0]);
    const T cos_value = cos(x.d[0]);
    result.d[0] = sin_value;
    for (int i = 1; i <= N; ++i)
      result.d[i] = x.d[i] * cos_value;
    for (int i = 0, k = 1 + N; i < N; ++i) {
      result.d[k] = (-.5 * sin_value) * sqr(x.d[i + 1]) + cos_value * x.d[k];
      ++k;
      for (int j = i + 1; j < N; ++j, ++k)
        result.d[k] = -sin_value * x.d[i + 1] * x.d[j + 1] + cos_value * x.d[k];
    }
    return result;
  }

  friend inline This cos(const This &x) {
    using std::sin;
    using std::cos;
    This result;
    const T sin_value = sin(x.d[0]);
    const T cos_value = cos(x.d[0]);
    result.d[0] = cos_value;
    for (int i = 1; i <= N; ++i)
      result.d[i] = x.d[i] * -sin_value;
    for (int i = 0, k = 1 + N; i < N; ++i) {
      result.d[k] = (-.5 * cos_value) * sqr(x.d[i + 1]) - sin_value * x.d[k];
      ++k;
      for (int j = i + 1; j < N; ++j, ++k)
        result.d[k] = -cos_value * x.d[i + 1] * x.d[j + 1] - sin_value * x.d[k];
    }
    return result;
  }

  friend inline This exp(const This &x) {
    using std::exp;
    This result;
    const T exp_value = exp(x.d[0]);
    result.d[0] = exp_value;
    for (int i = 1; i <= N; ++i)
      result.d[i] = x.d[i] * exp_value;
    for (int i = 0, k = 1 + N; i < N; ++i) {
      result.d[k] = exp_value * (.5 * sqr(x.d[i + 1]) + x.d[k]);
      ++k;
      for (int j = i + 1; j < N; ++j, ++k)
        result.d[k] = exp_value * (x.d[i + 1] * x.d[j + 1] + x.d[k]);
    }
    return result;
  }

  friend inline This log(const This &x) {
    using std::log;
    This result;
    result.d[0] = log(x.d[0]);
    const T inv = inverse(x.d[0]);
    for (int i = 1; i <= N; ++i)
      result.d[i] = inv * x.d[i];
    for (int i = 0, k = 1 + N; i < N; ++i) {
      result.d[k] = inv * (x.d[k] - .5 * inv * sqr(x.d[i + 1]));
      ++k;
      for (int j = i + 1; j < N; ++j, ++k)
        result.d[k] = inv * (x.d[k] - inv * x.d[i + 1] * x.d[j + 1]);
    }
    return result;
  }

  friend inline This coth(const This &x) {
    using std::sinh;
    This result;
    const T coth_value = coth(x.d[0]);
    const T eta = 1 - sqr(coth_value);
    // const T eta2 = -sqr(inverse(sinh(x.d[0])));  // This gets NaNs.
    result.d[0] = coth_value;
    for (int i = 1; i <= N; ++i)
      result.d[i] = x.d[i] * eta;
    for (int i = 0, k = 1 + N; i < N; ++i) {
      result.d[k] = eta * (x.d[k] - coth_value * sqr(x.d[i + 1]));
      ++k;
      for (int j = i + 1; j < N; ++j, ++k)
        result.d[k] = eta * (x.d[k] - 2 * coth_value * x.d[i + 1] * x.d[j + 1]);
    }
    return result;
  }

  friend inline bool isnan(const This &x) {
    for (int i = 0; i < M_; ++i)
      if (isnan(x.d[i]))
        return true;
    return false;
  }
};

#endif
