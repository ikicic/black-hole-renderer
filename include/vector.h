#ifndef VECTOR_H
#define VECTOR_H

#include <cmath>
#include <iostream>

#define _RESULT_FOR(i) \
  Vector result; \
  for (int i = 0; i < _N; ++i)

#define _FOR(i) for (int i = 0; i < _N; ++i)

template<typename _T, int _N> struct Vector {
  typedef _T value_type;

  _T v[_N];

  inline _T &operator[](int x) {
    return v[x];
  }
  inline const _T &operator[](int x) const {
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

  inline friend Vector operator*(_T x, const Vector &b) {
    _RESULT_FOR(i)
      result.v[i] = x * b.v[i];
    return result;
  }

  inline friend Vector operator*(const Vector &a, _T x) {
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

  inline Vector& operator*=(const _T &x) {
    _FOR(i)
      v[i] *= x;
    return *this;
  }

  inline Vector& operator/=(const _T &x) {
    _FOR(i)
      v[i] /= x;
    return *this;
  }

  inline Vector operator-() const {
    _RESULT_FOR(i)
      result.v[i] = -v[i];
    return result;
  }

  inline _T sqr_length(void) const {
    _T sum = _T();
    _FOR(i)
      sum += v[i] * v[i];
    return sum;
  }
  inline _T length(void) const {
    return sqrt(sqr_length());
  }
  inline Vector& normalize(void) {
    return *this /= length();
  }

  friend std::ostream& operator<<(std::ostream& stream, const Vector &A) {
    stream << '(';
    for (int i = 0; i < _N; ++i)
      stream << A.v[i] << (i == _N - 1 ? ')' : ' ');
    return stream;
  }

  template <typename _Vector>
  inline auto dot(const _Vector &a) const {
    typedef decltype(_T() * _Vector()[0]) result_type;
    result_type sum = result_type();
    _FOR(i)
      sum += v[i] * a.v[i];
    return sum;
  }
};


template <typename _T1, typename _T2>
inline Vector<decltype(_T1() * _T2()), 3> cross(const Vector<_T1, 3> &A,
                                                const Vector<_T2, 3> &B) {
  return {{
    A.v[1] * B.v[2] - A.v[2] * B.v[1],
    A.v[2] * B.v[0] - A.v[0] * B.v[2],
    A.v[0] * B.v[1] - A.v[1] * B.v[0]
  }};
}

template <typename _T, int _N>
Vector<_T, _N + 1> extend_vector(const _T &v0, const Vector<_T, _N> &V) {
  Vector<_T, _N + 1> result;
  result[0] = v0;
  for (int i = 0; i < _N; ++i)
    result[i + 1] = V[i];
  return result;
}


// utility.h
template <typename _T, int _N>
inline auto numerical_sqr_distance(
    const Vector<_T, _N> &A, const Vector<_T, _N> &B) {
  return (A - B).sqr_length();
}


#undef _RESULT_FOR
#undef _FOR
#endif
