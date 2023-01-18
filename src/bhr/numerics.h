#ifndef NUMERIC_H
#define NUMERIC_H

namespace bhr {

template <typename T, typename Coef>
T evaluate_polynomial(const Coef *coef, int N, const T &x) {
  T result = T(coef[N - 1]);
  for (int i = N - 2; i >= 0; --i)
    result = x * result + coef[i];
  return result;
}

template <typename T, typename Coef, int N>
inline T evaluate_polynomial(const Coef (&coef)[N], const T &x) {
  return evaluate_polynomial(coef, N, x);
}

template <typename T, typename Coef>
inline T evaluate_pade(
    const Coef *num, int N, const Coef *den, int M, const T &x) {
  // SPEED. These two polynomials can be evaluated simultaneously.
  return evaluate_polynomial(num, N, x) / evaluate_polynomial(den, M, x);
}

template <typename T, typename Coef, int N, int M>
inline T evaluate_pade(
    const Coef (&num)[N], const Coef (&den)[M], const T &x) {
  return evaluate_pade(num, N, den, M, x);
}

}  // namespace bhr

#endif
