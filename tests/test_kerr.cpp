#include "tests.h"
#if TESTS_ENABLED
#include "../include/coordinate.h"
#include "../include/kerr.h"
#include "../include/autodiff.h"

template <typename _T, typename _Coord>
void _get_kerr_autodiff_christoffel(
    const KerrSpacetime &kerr,
    const _Coord &position,
    Christoffel<_T> &output) {

  Matrix4<_T> invmetric;
  Matrix4<first_partial_derivatives<_T, 4> > metric_ad;
  BoyerLindquistVector4<first_partial_derivatives<_T, 4>> position_ad;
  for (int i = 0; i < 4; ++i) {
    position_ad[i].value() = position[i];
    for (int j = 0; j < 4; ++j)
      position_ad[i].first(j) = i == j;
  }
  metric_ad = kerr.get_metric_ll(position_ad);

  Matrix4<_T> metric;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      metric[i][j] = metric_ad[i][j].value();
  invmetric = matrix4_inverse(metric);


  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      for (int k = j; k < 4; ++k) {
        double sum = 0;
        for (int l = 0; l < 4; ++l) {
          sum += invmetric[i][l] * (
              metric_ad[l][j].first(k)
              + metric_ad[l][k].first(j)
              - metric_ad[j][k].first(l));
        }
        output[i][k][j] = output[i][j][k] = sum / 2;
      }
    }
  }
}

bool test_kerr(void) {
  Christoffel<double> christoffel1;
  Christoffel<double> christoffel2;
#if PREDEFINED_PARAMS
  KerrSpacetime kerr;
#else
  KerrSpacetime kerr(BLACK_HOLE_M, BLACK_HOLE_a);
#endif
  BoyerLindquistVector4<double> position{0, 5.123, 1.25, 1.53};

  _get_kerr_autodiff_christoffel(kerr, position, christoffel1);
  christoffel2 = kerr.get_christoffel_ull(position);

  BoyerLindquistVector4<double> direction, result1, result2, result3;
  for (int _test = 0; _test < 100; ++_test) {
    for (int k = 0; k < 4; ++k)
      direction[k] = random_double(-1., 1.);

    result1 = {0, 0, 0, 0};
    result2 = {0, 0, 0, 0};
    for (int k = 0; k < 4; ++k)
      for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b) {
          result1[k] += direction[a] * direction[b] * christoffel1[k][a][b];
          result2[k] += direction[a] * direction[b] * christoffel2[k][a][b];
        }
    result3 = kerr.geodesic_acceleration(position, direction);

    for (int k = 0; k < 4; ++k) {
      if (std::abs(result1[k] - result2[k]) > 1e-5
          || std::abs(result1[k] - result2[k]) > 1e-5) {
        std::cerr << "            Direction: " << direction << "\n";
        std::cerr << "             Autodiff: " << result1 << "\n";
        std::cerr << "  get_christoffel_ull: " << result2 << "\n";
        std::cerr << "geodesic_acceleration: " << result3 << "\n";
        for (int l = 0; l < 4; ++l) {
          for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j)
              fprintf(stderr, "%d%d%d %10.5lf  ", l, i, j, christoffel1[l][i][j] - christoffel2[l][i][j]);
            fprintf(stderr, "\n");
          }
        }

        // fprintf(stderr, "Christoffel [%d, %d, %d] Expected: %lg Received: %lg\n",
        //     k, i, j, expected, accel[k]);
        return false;
      }
    }
  }

  return true;
}
#endif
