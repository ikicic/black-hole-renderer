#include "base.h"
#include "parameters.h"
#include "qed_lagrangian.h"

void debug_units(void) {
#define PRINT(x) fprintf(stderr, "%20s = %12lg\n", #x, double(x));
  PRINT(UNIT_kg);
  PRINT(UNIT_m);
  PRINT(UNIT_s);
  PRINT(UNIT_C);
  PRINT(UNIT_T);
  PRINT(PHY_hbar);
  PRINT(BLACK_HOLE_M);
  PRINT(BLACK_HOLE_a);
  PRINT(BLACK_HOLE_Q);
  PRINT(_BLACK_HOLE_a_MAX);
  PRINT(_BLACK_HOLE_Q_MAX);
  PRINT(_BLACK_HOLE_r_S);
  fprintf(stderr, "%20s = %12lg s^2 C^4 kg^-3 m^-1\n", "QED::_b",
      QED::_b / (sqr(UNIT_s * sqr(UNIT_C)) / cube(UNIT_kg) / UNIT_m));
  fprintf(stderr, "%20s = %12lg T^-2\n", "lambda1", QED::lambda1 * sqr(UNIT_T));
  fprintf(stderr, "%20s = %12lg T^-2\n", "lambda2", QED::lambda2 * sqr(UNIT_T));
#if RENDER_DISK
  PRINT(_KERR_ISCO);
  PRINT(INNER_RADIUS);
  PRINT(OUTER_RADIUS);
#endif
  PRINT(PHY_Bc);
  fprintf(stderr, "%20s = %12lg km\n", "Black hole r_S",
      _BLACK_HOLE_r_S / UNIT_km);
  fprintf(stderr, "%20s = %12lg km = %12lg\n", "Neutron star r",
      NEUTRON_STAR_r / UNIT_km, NEUTRON_STAR_r);
  // assert(false);
#undef PRINT
}
