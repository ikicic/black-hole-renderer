#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <bhr/base.h>
#include <bhr/physical_constants.h>

namespace bhr {

static constexpr QUANTITY(0, 1, 0, 0) _BLACK_HOLE_r_S = 2 * PHY_G * BLACK_HOLE_M / sqr(PHY_c);
static constexpr QUANTITY(0, 1, 0, 0) _BLACK_HOLE_a_MAX = _BLACK_HOLE_r_S / 2;
static constexpr QUANTITY(0, 0, 0, 1) _BLACK_HOLE_Q_MAX =
    constexpr_sqrt(PHY_G / PHY_Coulomb_k) * BLACK_HOLE_M;
static constexpr QUANTITY(0, 1, 0, 0) BLACK_HOLE_a = 0.998 * _BLACK_HOLE_a_MAX;
static constexpr QUANTITY(0, 0, 0, 1) BLACK_HOLE_Q = 0.1 * _BLACK_HOLE_Q_MAX;

static CMATH_CONSTEXPR QUANTITY(0, 1, 0, 0) _BLACK_HOLE_radius =
    BLACK_HOLE_M + constexpr_sqrt(sqr(BLACK_HOLE_M) - sqr(BLACK_HOLE_a));

static constexpr QUANTITY(0, 1, 0, 0) MAX_r = _BLACK_HOLE_r_S * 10000;

template <typename T>
CMATH_CONSTEXPR auto kerr_direct_ISCO(T M, T a) {
  auto rs = (2 * PHY_G / sqr(PHY_c)) * M;
  auto chi = 2 * a / rs;
  auto t1 = std::cbrt(1 - chi);
  auto t2 = std::cbrt(1 + chi);
  auto Z1 = 1 + t1 * t2 * (t1 + t2);
  auto Z2 = std::sqrt(3 * sqr(chi) + sqr(Z1));
  return (rs / 2) * (3 + Z2 - std::sqrt((3 - Z1) * (2 + Z1 + 2 * Z2)));
}

#if RENDER_DISK
static CMATH_CONSTEXPR QUANTITY(0, 1, 0, 0) _KERR_ISCO = kerr_direct_ISCO(BLACK_HOLE_M, BLACK_HOLE_a);
#endif

#if RENDER_DISK == DISK_KERTAP
static CMATH_CONSTEXPR QUANTITY(0, 1, 0, 0) INNER_RADIUS = _KERR_ISCO;
static constexpr QUANTITY(0, 1, 0, 0) OUTER_RADIUS = 10 * _BLACK_HOLE_r_S;
#elif RENDER_DISK == DISK_DUMMY
static constexpr QUANTITY(0, 1, 0, 0) INNER_RADIUS = 3 * _BLACK_HOLE_r_S;
static constexpr QUANTITY(0, 1, 0, 0) OUTER_RADIUS = 10 * _BLACK_HOLE_r_S;
#elif RENDER_DISK == DISK_SHAKURA
static CMATH_CONSTEXPR QUANTITY(0, 1, 0, 0) INNER_RADIUS = _KERR_ISCO;
static constexpr QUANTITY(0, 1, 0, 0) OUTER_RADIUS = 100 * _BLACK_HOLE_r_S;
#endif

// QUANTITY(0, 1, 0, 0) INNER_RADIUS = 6 * _BLACK_HOLE_r_S;
// QUANTITY(0, 1, 0, 0) OUTER_RADIUS = 10 * INNER_RADIUS;

}  // namespace bhr

#endif
