#ifndef NUMERIC_H
#define NUMERIC_H


template <typename _T, typename _Coef>
_T evaluate_polynomial(const _Coef *coef, int N, const _T &x) {
  _T result = _T(coef[N - 1]);
  for (int i = N - 2; i >= 0; --i)
    result = x * result + coef[i];
  return result;
}

template <typename _T, typename _Coef, int N>
inline _T evaluate_polynomial(const _Coef (&coef)[N], const _T &x) {
  return evaluate_polynomial(coef, N, x);
}

template <typename _T, typename _Coef>
inline _T evaluate_pade(
    const _Coef *num, int N, const _Coef *den, int M, const _T &x) {
  // SPEED. These two polynomials can be evaluated simultaneously.
  return evaluate_polynomial(num, N, x) / evaluate_polynomial(den, M, x);
}

template <typename _T, typename _Coef, int N, int M>
inline _T evaluate_pade(
    const _Coef (&num)[N], const _Coef (&den)[M], const _T &x) {
  return evaluate_pade(num, N, den, M, x);
}

#endif
