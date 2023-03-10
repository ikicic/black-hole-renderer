#ifndef PHYSICAL_CONSTANTS_H
#define PHYSICAL_CONSTANTS_H

#include <cassert>
#include <bhr/float_helpers.h>
#include <bhr/utility.h>

namespace bhr {

#define TEST_UNITS  0
#if TEST_UNITS
# define QUANTITY(a, b, c, d) Quantity<(a), (b), (c), (d)>
#else
# define QUANTITY(a, b, c, d) double
#endif


/* Mass, Length, Time, Charge */
template <int A, int B, int C, int D>
struct Quantity {
  const double value;

  constexpr explicit Quantity(double _value) : value(_value) {
    // static_assert(A == 0 && B == 0 && C == 0 && D == 0);
  }

  friend inline constexpr Quantity operator+(Quantity x, Quantity y) {
    return Quantity(x.value + y.value);
  }

  friend inline constexpr Quantity operator-(Quantity x, Quantity y) {
    return Quantity(x.value - y.value);
  }
  friend inline constexpr Quantity operator*(double x, Quantity y) {
    return Quantity(x * y.value);
  }
  friend inline constexpr Quantity operator*(Quantity x, double y) {
    return Quantity(x.value * y);
  }
  friend inline constexpr Quantity operator/(Quantity x, double y) {
    return Quantity(x.value / y);
  }
  friend inline constexpr auto operator/(double x, Quantity y) {
    return Quantity<-A, -B, -C, -D>(x / y.value);
  }

  explicit inline constexpr operator double(void) const {
    return value;
  }
};

template <int A, int B, int C, int D>
inline constexpr auto sqrt(Quantity<A, B, C, D> x) {
  using std::sqrt;
  static_assert(A % 2 == 0 && B % 2 == 0 && C % 2 == 0 && D % 2 == 0);
  return Quantity<A / 2, B / 2, C / 2, D / 2>(sqrt(x.value));
}

template <int A, int B, int C, int D>
inline constexpr auto constexpr_sqrt(Quantity<A, B, C, D> x) {
  static_assert(A % 2 == 0 && B % 2 == 0 && C % 2 == 0 && D % 2 == 0);
  return Quantity<A / 2, B / 2, C / 2, D / 2>(constexpr_sqrt(x.value));
}


template <int A1, int B1, int C1, int D1,
          int A2, int B2, int C2, int D2>
inline constexpr auto operator*(
    Quantity<A1, B1, C1, D1> x, Quantity<A2, B2, C2, D2> y) {
  return Quantity<A1 + A2, B1 + B2, C1 + C2, D1 + D2>(
      x.value * y.value);
}

template <int A1, int B1, int C1, int D1,
          int A2, int B2, int C2, int D2>
inline constexpr auto operator/(
    Quantity<A1, B1, C1, D1> x, Quantity<A2, B2, C2, D2> y) {
  return Quantity<A1 - A2, B1 - B2, C1 - C2, D1 - D2>(
      x.value / y.value);
}

//       M,  L,  T, Q
constexpr QUANTITY(-1, 3, -2, 0) PHY_G{1};
constexpr QUANTITY(0, 1, -1, 0) PHY_c{1};
constexpr QUANTITY(1, 1, 0, -2) PHY_mu0{1};
constexpr QUANTITY(-1, -3, 2, 2) PHY_eps0 = 1. / (PHY_mu0 * PHY_c * PHY_c);
constexpr QUANTITY(1, 3, -2, -2) PHY_Coulomb_k = 1. / (4 * M_PI * PHY_eps0);

constexpr QUANTITY(1, 0, 0, 0) BLACK_HOLE_M{0.05};

// http://adsabs.harvard.edu/abs/2013HEAD...1310308D
constexpr  QUANTITY(1, 0, 0, 0) PHY_M_Sun = BLACK_HOLE_M / 1.25;
// constexpr  QUANTITY(1, 0, 0, 0) PHY_M_Sun = BLACK_HOLE_M / 10;

constexpr QUANTITY(1, 0, 0, 0) UNIT_kg = PHY_M_Sun / 1.989e30;
constexpr QUANTITY(0, 1, 0, 0) UNIT_m =
    UNIT_kg * (PHY_G / 6.67408e-11) / sqr(PHY_c / 299792458.);
constexpr QUANTITY(0, 1, 0, 0) UNIT_km = 1000. * UNIT_m;

// constexpr QUANTITY(0, 1, 0, 0) NEUTRON_STAR_r = 12.4 * UNIT_km;
constexpr QUANTITY(0, 1, 0, 0) NEUTRON_STAR_r = 0 * UNIT_km;

constexpr QUANTITY(0, 0, 1, 0) UNIT_s = 299792458. * UNIT_m / PHY_c;
constexpr QUANTITY(0, 0, 1, 0) UNIT_y = 365.24 * 86400 * UNIT_s;
constexpr QUANTITY(0, 0, -1, 0) UNIT_Hz = 1 / UNIT_s;
constexpr QUANTITY(1, 1,-2, 0) UNIT_N = UNIT_kg * UNIT_m / UNIT_s / UNIT_s;
constexpr QUANTITY(1, 2,-2, 0) UNIT_J = UNIT_N * UNIT_m;
constexpr QUANTITY(1, 2,-3, 0) UNIT_W = UNIT_J / UNIT_s;
constexpr QUANTITY(0, 0,-2, 2) UNIT_A_sqr = 4e-7 * M_PI * UNIT_N / PHY_mu0;
constexpr QUANTITY(0, 0,-1, 1) UNIT_A = constexpr_sqrt(UNIT_A_sqr);
constexpr QUANTITY(0, 0, 0, 1) UNIT_C = UNIT_A * UNIT_s;

constexpr QUANTITY(1, 2,-2,-1) UNIT_V = UNIT_kg * sqr(UNIT_m) / (sqr(UNIT_s) * UNIT_C);
constexpr QUANTITY(1, 0,-1,-1) UNIT_T = UNIT_kg / (UNIT_A * UNIT_s * UNIT_s);
constexpr QUANTITY(1, 0,-1,-1) UNIT_G = 1e-4 * UNIT_T;

constexpr QUANTITY(1, 2,-1, 0) PHY_h = 6.62607004e-34 * UNIT_m * UNIT_m * UNIT_kg / UNIT_s;
constexpr QUANTITY(1, 2,-1, 0) PHY_hbar = PHY_h / (2 * M_PI);
constexpr QUANTITY(0, 0, 0, 1) PHY_e = 1.602176565e-19 * UNIT_C;
constexpr QUANTITY(1, 0, 0, 0) PHY_m_e = 9.10938356e-31 * UNIT_kg;

constexpr QUANTITY(0, 0, 0, 0) PHY_alpha(7.2973525664e-3);
constexpr QUANTITY(1, 0,-1,-1) PHY_Bc = sqr(PHY_m_e * PHY_c) / PHY_e / PHY_hbar;
constexpr QUANTITY(2, 0,-2,-2) PHY_sqr_Bc = sqr(PHY_Bc);
constexpr QUANTITY(-2, 0, 2, 2) PHY_inv_sqr_Bc = 1 / PHY_sqr_Bc;


constexpr QUANTITY(0, 0, 0, 0) UNIT_K = 1.0;  // NO DIM CHECK!
constexpr QUANTITY(1, 2,-2, 0) PHY_kB = 1.38064852e-23 * UNIT_J / UNIT_K;
constexpr QUANTITY(1, 2,-2, 0) PHY_black_body_sigma =
    (2 * M_PI * sqr(sqr(M_PI)) / 15) * sqr(sqr(PHY_kB) / PHY_c) / cube(PHY_h);

void debug_units(void);

}  // namespace bhr

#endif
