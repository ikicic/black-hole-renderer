#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <bhr/utility.h>


// https://math.okstate.edu/people/yqwang/teaching/math4513_fall11/Notes/rungekutta.pdf
// Adaptive Runge-Kutta-FehlBerg
template <typename _Vector, typename _RHS>
inline real_t integration_step__RGF45(
    const _RHS &RHS, const real_t h, const _Vector &u, _Vector *out) {

#define EVALUATE(v, w) { \
      w = RHS((v)); \
      mult(&(w), h); \
    }
#define FRAC(a, b) (real_t(a) / b)
  _Vector k1, k2, k3, k4, k5, k6, v;

  EVALUATE(u, k1);
  set_and_mult_add(&v, u, FRAC(1, 4), k1);
  EVALUATE(v, k2);
  set_and_mult_add(&v, u, FRAC(3, 32), k1, FRAC(9, 32), k2);
  EVALUATE(v, k3);
  set_and_mult_add(&v, u,
      FRAC(1932, 2197), k1,
      FRAC(-7200, 2197), k2,
      FRAC(7296, 2197), k3);
  EVALUATE(v, k4);
  set_and_mult_add(&v, u,
      FRAC(439, 216), k1,
      FRAC(-8, 1), k2,
      FRAC(3680, 513), k3,
      FRAC(-845, 4104), k4);
  EVALUATE(v, k5);
  set_and_mult_add(&v, u,
      FRAC(-8, 27), k1,
      FRAC(2, 1), k2,
      FRAC(-3544, 2565), k3,
      FRAC(1859, 4104), k4,
      FRAC(-11, 40), k5);
  EVALUATE(v, k6);


  set_and_mult_add(&v, u,
      FRAC(25, 216), k1,
      FRAC(1408, 2565), k3,
      FRAC(2197, 4104), k4,
      FRAC(-1, 5), k5);
  set_and_mult_add(out, u,
      FRAC(16, 135), k1,
      FRAC(6656, 12825), k3,
      FRAC(28561, 56430), k4,
      FRAC(-9, 50), k5,
      FRAC(2, 55), k6);

#undef EVALUATE
#undef FRAC

  // Not sure which norm should be used here...
  // R = |v - *u| / h
  return numerical_distance(v, *out) / h;
}


template <typename _Vector, typename _RHS>
inline void integrate__RGF45(
    const _RHS &RHS,
    _Vector initial,
    real_t t,
    const real_t min_h,
    real_t h,  // start_h
    const real_t max_h,
    const real_t epsilon,
    _Vector * const output) {

  constexpr real_t SAFETY = real_t(0.84);
  _Vector out;

  while (t != 0.0) {
    if (h > t) h = t;
    const real_t R = integration_step__RGF45(RHS, h, initial, &out);
    if (R <= epsilon) {
      t -= h;
      initial = out;
    } else if (h <= min_h) {
      t -= h;
      initial = out;
      continue;
    }
    const real_t delta = SAFETY * std::pow(epsilon / R, 0.25);
    h = clamp(h * delta, min_h, max_h);
  }

  *output = out;
}
#endif
