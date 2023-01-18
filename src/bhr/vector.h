#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>
#include <iostream>

#define _RESULT_FOR(i) \
  Vector result; \
  for (int i = 0; i < N; ++i)

#define _FOR(i) for (int i = 0; i < N; ++i)

template<typename T, int N> struct Vector {
  typedef T value_type;

  T v[N];

  inline T &operator[](int x) {
    return v[x];
  }
  inline const T &operator[](int x) const {
    return v[x];
  }

  inline friend Vector operator+(const Vector &a, const Vector &b) {
    _RESULT_FOR(i)
      result.v[i] = a.v[i] + b.v[i];
    return result;
  }

  inline friend Vector operator-(const Vector &a, const Vector &b) {
    _RESULT_FOR(i)
      result.v[i] = a.v[i] - b.v[i];
    return result;
  }

  inline friend Vector operator*(T x, const Vector &b) {
    _RESULT_FOR(i)
      result.v[i] = x * b.v[i];
    return result;
  }

  inline friend Vector operator*(const Vector &a, T x) {
    _RESULT_FOR(i)
      result.v[i] = a.v[i] * x;
    return result;
  }

  inline Vector& operator+=(const Vector &a) {
    _FOR(i)
      v[i] += a.v[i];
    return *this;
  }

  inline Vector& operator-=(const Vector &a) {
    _FOR(i)
      v[i] -= a.v[i];
    return *this;
  }

  inline Vector& operator*=(const T &x) {
    _FOR(i)
      v[i] *= x;
    return *this;
  }

  inline Vector& operator/=(const T &x) {
    _FOR(i)
      v[i] /= x;
    return *this;
  }

  inline Vector operator-() const {
    _RESULT_FOR(i)
      result.v[i] = -v[i];
    return result;
  }

  inline T sqr_length(void) const {
    T sum = T();
    _FOR(i)
      sum += v[i] * v[i];
    return sum;
  }
  inline T length(void) const {
    return sqrt(sqr_length());
  }
  inline Vector& normalize(void) {
    return *this /= length();
  }

  friend std::ostream& operator<<(std::ostream& stream, const Vector &A) {
    stream << '(';
    for (int i = 0; i < N; ++i)
      stream << A.v[i] << (i == N - 1 ? ')' : ' ');
    return stream;
  }

  template <typename Vector>
  inline auto dot(const Vector &a) const {
    typedef decltype(T() * Vector()[0]) result_type;
    result_type sum = result_type();
    _FOR(i)
      sum += v[i] * a.v[i];
    return sum;
  }
};


template <typename T1, typename T2>
inline Vector<decltype(T1() * T2()), 3> cross(const Vector<T1, 3> &A,
                                                const Vector<T2, 3> &B) {
  return {{
    A.v[1] * B.v[2] - A.v[2] * B.v[1],
    A.v[2] * B.v[0] - A.v[0] * B.v[2],
    A.v[0] * B.v[1] - A.v[1] * B.v[0]
  }};
}

template <typename T, int N>
Vector<T, N + 1> extend_vector(const T &v0, const Vector<T, N> &V) {
  Vector<T, N + 1> result;
  result[0] = v0;
  for (int i = 0; i < N; ++i)
    result[i + 1] = V[i];
  return result;
}


// utility.h
template <typename T, int N>
inline auto numerical_sqr_distance(
    const Vector<T, N> &A, const Vector<T, N> &B) {
  return (A - B).sqr_length();
}


#undef _RESULT_FOR
#undef _FOR
#endif
