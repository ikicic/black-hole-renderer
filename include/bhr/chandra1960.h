#ifndef CHANDRA_1960
#define CHANDRA_1960

#include <utility>

inline std::pair<double, double> chandra1960(double mu) {
  // Copied from KERTAP::chandra60.m.

  // delta is the degree of polarization
  // weight is the limbdarkening factor.
  // the table used of interpolation comes from Chandra1960 Tab XXIV

  static const double I_list[] = {
      0.41441, 0.47490, 0.52397, 0.57001, 0.61439, 0.65770,
      0.70029, 0.74234, 0.78398, 0.82530, 0.86637, 0.90722,
      0.94789, 0.98842, 1.02882, 1.06911, 1.10931, 1.14943,
      1.18947, 1.22945, 1.26938
  };
  static const double del_list[] = {
      0.11713, 0.08979, 0.07448, 0.06311, 0.05410, 0.04667,
      0.04041, 0.03502, 0.03033, 0.02619, 0.02252, 0.01923,
      0.01627, 0.01358, 0.011123, 0.008880, 0.006818, 0.004919,
      0.003155, 0.001522, 0
  };
  constexpr int N = sizeof(I_list) / sizeof(I_list[0]);
  constexpr int M = sizeof(del_list) / sizeof(del_list[0]);

  return std::make_pair(
      linear_interpolation(I_list, N, 0., mu, 1.),
      linear_interpolation(del_list, M, 0., mu, 1.)
  );
}

#endif
