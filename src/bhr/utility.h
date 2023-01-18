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

template<typename T> constexpr inline T inverse(const T &x) {
  /* How do you properly make the multiplicative inverse? */
  return T(1) / x;
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

template <typename Out, typename T>
inline Out linear_interpolation(
    const Out *input, int N, const T &low, const T &x, const T &high) {
  if (x <= low)
    return input[0];
  T scaled = (N - 1) * (x - low) / (high - low);
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

template <typename T>
inline auto numerical_sqr_distance(const T &A, const T &B) {
  return sqr(A - B);
}

template <typename T>
inline auto numerical_distance(const T &A, const T &B) {
  return std::sqrt(numerical_sqr_distance(A, B));
}

template <typename T>
inline auto numerical_sqr_magnitude(const T &A) {
  return sqr(A);
}

template <typename T>
inline T clamp(const T &value, const T &low, const T &high) {
  return value > low ? (value < high ? value : high) : low;
}

inline bool ends_with(const std::string &s, const std::string &suffix) {
  return s.size() >= suffix.size() &&
      s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}




/* Optimized arithmetics operations */

/* Default operations. */
template <typename State, typename T>
inline void mult(State *lhs, const T &factor) {
  *lhs *= factor;
}

template <typename State, typename T>
inline void mult_add(State *lhs, const T &factor, const State &A) {
  *lhs += factor * A;
}

template <typename State, typename T>
inline void set_and_mult_add(
    State *lhs, const State &A, const T &factor, const State &B) {
  *lhs = A + factor * B;
}

/* Default multi operations. */
template <typename Vector>
inline void mult_add_multiple(Vector * /*out*/) {
  // No-op.
}

template <typename Vector, typename T, typename ...Args>
inline void mult_add_multiple(
    Vector *out, const T &c, const Vector &A, Args ...args) {
  mult_add(out, c, A);
  mult_add_multiple(out, args...);
}

template <typename State, typename T, typename ...Args>
inline void set_and_mult_add(State *out, const State &A, const T &c, const State &B, Args ...args) {
  set_and_mult_add(out, A, c, B);
  mult_add_multiple(out, args...);
}


#endif
