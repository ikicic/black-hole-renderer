#ifndef UTILITY_H
#define UTILITY_H

#include <cassert>

#ifndef __clang__
# define CMATH_CONSTEXPR  constexpr
#else
# define CMATH_CONSTEXPR
#endif

template <typename T>
constexpr inline auto sqr(T x) { return x * x; }

template <typename T>
constexpr inline auto cube(T x) { return sqr(x) * x; }

template<typename _T> constexpr inline _T inverse(const _T &x) {
  /* How do you properly make the multiplicative inverse? */
  return _T(1) / x;
}

template <typename T>
constexpr inline auto inv_sqrt(const T &x) {
  using std::sqrt;
  return inverse(sqrt(x));
}

template <typename T>
constexpr inline auto coth(T x) {
  using std::tanh;
  return inverse(tanh(x));
}

template <typename _Out, typename _T>
inline _Out linear_interpolation(
    const _Out *input, int N, const _T &low, const _T &x, const _T &high) {
  if (x <= low)
    return input[0];
  _T scaled = (N - 1) * (x - low) / (high - low);
  int index = int(scaled);
  // assert(0 <= index);
  // assert(index < N);
  if (index >= N - 1)
    return input[N - 1];
  return (1 - (scaled - index)) * input[index]
       + (scaled - index) * input[index + 1];
}


inline double numerical_distance(double A, double B) {
  return std::fabs(A - B);
}

template <typename _T>
inline auto numerical_sqr_distance(const _T &A, const _T &B) {
  return sqr(A - B);
}

template <typename _T>
inline auto numerical_distance(const _T &A, const _T &B) {
  return std::sqrt(numerical_sqr_distance(A, B));
}

template <typename _T>
inline auto numerical_sqr_magnitude(const _T &A) {
  return sqr(A);
}

template <typename _T>
inline _T clamp(const _T &value, const _T &low, const _T &high) {
  return value > low ? (value < high ? value : high) : low;
}

inline bool ends_with(const std::string &s, const std::string &suffix) {
  return s.size() >= suffix.size() &&
      s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}




/* Optimized arithmetics operations */

/* Default operations. */
template <typename _State, typename _T>
inline void mult(_State *lhs, const _T &factor) {
  *lhs *= factor;
}

template <typename _State, typename _T>
inline void mult_add(_State *lhs, const _T &factor, const _State &A) {
  *lhs += factor * A;
}

template <typename _State, typename _T>
inline void set_and_mult_add(
    _State *lhs, const _State &A, const _T &factor, const _State &B) {
  *lhs = A + factor * B;
}

/* Default multi operations. */
template <typename _Vector>
inline void mult_add_multiple(_Vector * /*out*/) {
  // No-op.
}

template <typename _Vector, typename _T, typename ...Args>
inline void mult_add_multiple(
    _Vector *out, const _T &c, const _Vector &A, Args ...args) {
  mult_add(out, c, A);
  mult_add_multiple(out, args...);
}

template <typename _State, typename _T, typename ...Args>
inline void set_and_mult_add(_State *out, const _State &A, const _T &c, const _State &B, Args ...args) {
  set_and_mult_add(out, A, c, B);
  mult_add_multiple(out, args...);
}


#endif
