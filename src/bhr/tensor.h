#ifndef TENSOR_H
#define TENSOR_H

#include <bhr/autodiff_2nd.h>
#include <bhr/physical_constants.h>
#include <bhr/qed_lagrangian.h>

namespace bhr {

template <typename T> struct Christoffel {
  typedef T matrix44[4][4];
  matrix44 mat[4];

  inline matrix44 &operator[](int y) {
    return mat[y];
  }
  inline const matrix44 &operator[](int y) const {
    return mat[y];
  }

  friend std::ostream& operator<<(std::ostream& stream, const Christoffel &C) {
    stream << "[m^-1]\n";
    for (int k = 0; k < 4; ++k) {
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
          stream << C.mat[k][i][j] * UNIT_m << (j < 4 - 1 ? ' ' : '\n');
      stream << '\n';
    }
    return stream << '\n';
  }
};

inline void __print_matrix44(const double tmp[4][4]) {
  // fprintf(stderr, "--------------------------------------------------------\n");
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j)
      fprintf(stderr, "%12lg ", tmp[i][j]);
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
}

inline void __print_matrix444(const double tmp[4][4][4]) {
  // fprintf(stderr, "--------------------------------------------------------\n");
  for (int k = 0; k < 4; ++k) {
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j)
        fprintf(stderr, "%12lg ", tmp[k][i][j]);
      fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
}

template <
    typename T,
    template <typename> class Vector,
    typename MetricLLFunc,
    typename MetricUUFunc,
    typename PotentialLFunc>
inline Vector<T> geodesic_acceleration__magnetic_field__lowest_order(
    const MetricLLFunc &metric_ll_func,
    const MetricUUFunc &metric_uu_func,
    const PotentialLFunc &potential_l_func,
    const Vector<T> &position_u,
    const Vector<T> &vec_u) {

  constexpr double XI = QED::lambda2;

  typedef first_partial_derivatives<T, 4> fpds;
  const Vector<fpds> ad_position_u{
    fpds(position_u[0], 1, 0, 0, 0),
    fpds(position_u[1], 0, 1, 0, 0),
    fpds(position_u[2], 0, 0, 1, 0),
    fpds(position_u[3], 0, 0, 0, 1)
  };

  const Matrix4<fpds> metric_ll = metric_ll_func(ad_position_u);
  Christoffel<T> chr_lll; // Without the factor 1/2.
  for (int k = 0; k < 4; ++k)
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) {
        chr_lll[k][i][j] = metric_ll[k][j].first(i)
                         + metric_ll[k][i].first(j)
                         - metric_ll[i][j].first(k);
      }

  typedef second_partial_derivatives<T, 4> spds;
  const Vector<spds> add_position_u{
    spds(position_u[0], 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    spds(position_u[1], 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    spds(position_u[2], 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    spds(position_u[3], 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  };
  const Vector<spds> A_l = potential_l_func(add_position_u);
  /* Magnetic field tensor perturbation is of the form:
   * g~^ab = g^ab + xi gB^ab
   * g~_ab = g_ab - xi gB_ab
   * gB^ab = F^ac F_c^b = -g^cd F^ac F^bd
   * gB_ab = F_ac F^c_b = -g_cd F_ac F_bd
   */
  Matrix4<T> F_ll;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      F_ll[i][j] = A_l[j].first(i) - A_l[i].first(j);

  /* F_der_lll[k][i][j] = d(F_ij) / dx^k */
  T F_der_lll[4][4][4];
  for (int k = 0; k < 4; ++k)
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) {
        F_der_lll[k][i][j] = A_l[j].second(i, k)
                           - A_l[i].second(j, k);
      }

  /* gB_der[k][i][j] = d(g_ij) / dx^k */
  T gB_der_lll[4][4][4];
  for (int k = 0; k < 4; ++k)
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) {
        T tmp = T();
        for (int a = 0; a < 4; ++a)
          for (int b = 0; b < 4; ++b) {
            tmp += metric_ll[a][b].first(k) * F_ll[i][a] * F_ll[j][b];
            tmp += metric_ll[a][b].value() * (
                  F_der_lll[k][i][a] * F_ll[j][b]
                + F_der_lll[k][j][b] * F_ll[i][a]
            );
          }
        gB_der_lll[k][i][j] = -XI * tmp;  // gB_uu has + sign.
      }

  Christoffel<T> chr2_lll;
  for (int k = 0; k < 4; ++k)
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) {
        chr2_lll[k][i][j] = -(gB_der_lll[i][k][j]
                             + gB_der_lll[j][k][i]
                             - gB_der_lll[k][i][j]);
      }

  for (int k = 0; k < 4; ++k)
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        chr_lll[k][i][j] += chr2_lll[k][i][j];


  /* Calculate d(vec_k) / dlambda. */
  Vector<T> result_l = Vector<T>();
  for (int k = 0; k < 4; ++k)
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        result_l[k] += vec_u[i] * chr_lll[k][i][j] * vec_u[j];

  /* Raise index with the perturbed metric. */
  /* First, raise F_ll to F_uu and F_lu. */
  Matrix4<T> metric_uu = metric_uu_func(position_u);
  Matrix4<T> F_lu;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) {
      T tmp = T();
      for (int k = 0; k < 4; ++k)
        tmp += metric_uu[j][k] * F_ll[i][k];
      F_lu[i][j] = tmp;
    }

  Matrix4<T> F_uu__minus;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) {
      T tmp = T();
      for (int k = 0; k < 4; ++k) {
        // tmp += metric_uu[i][k] * (-F_lu[k][j]);
        tmp += metric_uu[i][k] * F_lu[j][k];
      }
      F_uu__minus[i][j] = tmp;
    }

  Matrix4<T> gB_uu;  // Metric magnetic field perturbation.
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) {
      T tmp = T();
      for (int k = 0; k < 4; ++k) {
        // tmp += F_uu[i][k] * F_lu[k][j];
        tmp += F_uu__minus[i][k] * F_lu[j][k];
      }
      gB_uu[i][j] = XI * tmp;
    }

  Vector<T> result_u;
  for (int i = 0; i < 4; ++i) {
    T tmp = T();
    for (int j = 0; j < 4; ++j)
      tmp += metric_uu[i][j] * result_l[j];
    T tmp2 = T();
    for (int j = 0; j < 4; ++j)
      tmp2 += gB_uu[i][j] * result_l[j];
    /* Christoffel symbols didn't include factor 1/2, so we include it here. */
    // fprintf(stderr, "%lg %lg\n", tmp, tmp2);
    result_u[i] = (tmp + tmp2) / 2;
  }

  return result_u;
}


template <
    typename T,
    template <typename> class Vector,
    typename MetricLLFunc,
    typename MetricUUFunc,
    typename FLLFunc,
    typename LagrangianFunc>
inline Vector<T> geodesic_acceleration__magnetic_field(
    const MetricLLFunc &metric_ll_func,
    const MetricUUFunc &metric_uu_func,
    const FLLFunc &F_ll_func,
    const LagrangianFunc &lagrangian_func,
    const Vector<T> &position_u,
    const Vector<T> &vec_u) {
  using std::sqrt;

  typedef first_partial_derivatives<T, 4> fpd4;
  const Vector<fpd4> ad_position_u{
    fpd4(position_u[0], 1, 0, 0, 0),
    fpd4(position_u[1], 0, 1, 0, 0),
    fpd4(position_u[2], 0, 0, 1, 0),
    fpd4(position_u[3], 0, 0, 0, 1)
  };
  Matrix4<fpd4> total_metric_ll = metric_ll_func(ad_position_u);
  Matrix4<fpd4> total_metric_uu = metric_uu_func(ad_position_u);

  const Matrix4<fpd4> F_ll = F_ll_func(ad_position_u);
  const Matrix4<fpd4> F_ul = total_metric_uu * F_ll;
  const Matrix4<fpd4> F_uu = F_ul * total_metric_uu;

  fpd4 F = fpd4();
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      F += F_ll[i][j] * F_uu[i][j];
  F /= 4;

  // // TODO: Check the factor.
  // fpd4 det = determinant4(total_metric_ll);
  // fpd4 G = (F_uu[0][1] * F_uu[2][3]
  //         - F_uu[0][2] * F_uu[1][3]
  //         + F_uu[0][3] * F_uu[1][2]) * sqrt(-det);
  fpd4 G(0, 0, 0, 0, 0);

  if (std::isnan(F.value()) || std::isnan(G.value())) {
    fpd4 det = determinant4(total_metric_ll);
    std::cerr << "NAN!!!\n";
    std::cerr << F << "  " << G << '\n';
    std::cerr << position_u << '\n';
    std::cerr << vec_u << '\n';
    std::cerr << det << '\n';
    exit(1);
  }
  if (debug >= 3) {
    std::cerr << "--------------------------------------------------------\n";
    std::cerr << "position=" << position_u << "  \t";
    std::cerr << "direction=" << vec_u << "  \t";
    std::cerr << "F=" << F / sqr(UNIT_T) << " T^2\t";
    std::cerr << "B=" << sqrt(2 * F) / UNIT_T << " T\n";
  }

  if (F + G > sqr(1e4 * UNIT_T)) {
      /* Calculate lambdas and d/dF and d/dG. */
    typedef first_partial_derivatives<T, 2> fpd2;
    fpd2 lambda2[2];
#if PREPROCESS_LAMBDAS
    using namespace QEDCache;
    (void)lagrangian_func;
    LambdaCacheType _cache = get_cached_lambdas(F.value(), G.value());
    lambda2[0] = _cache.lambda[0];
    lambda2[1] = _cache.lambda[1];
#else
    std::tie(lambda2[0], lambda2[1]) = qed_metric_correction_lambda(
        lagrangian_func, fpd2(F.value(), 1, 0), fpd2(G.value(), 0, 1));
#endif

    /* Convert d/dF and d/dG to d/dx^mu. */
    fpd4 lambda4[2];
    for (int i = 0; i < 2; ++i) {
      lambda4[i].value() = lambda2[i].value();
      for (int j = 0; j < 4; ++j) {
        lambda4[i].first(j) = lambda2[i].first(0) * F.first(j)
                            + lambda2[i].first(1) * G.first(j);
      }
    }

    total_metric_uu += lambda4[1] * (F_ul * F_uu);
    total_metric_ll = matrix4_inverse(total_metric_uu);

    if (debug >= 3) {
      std::cerr << "lambda[1] (at B=0)=" << QED::lambda2 << '\n';
      std::cerr << "lambda2[1]=" << lambda2[1] << '\n';
      std::cerr << "lambda4[1]=" << lambda4[1] << '\n';
    }
  }

  Christoffel<T> chr_lll; // Without the factor 1/2.
  for (int k = 0; k < 4; ++k)
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j) {
        chr_lll[k][i][j] = total_metric_ll[k][j].first(i)
                         + total_metric_ll[k][i].first(j)
                         - total_metric_ll[i][j].first(k);
      }
  if (debug >= 3) {
    std::cerr << "det=" << determinant4(total_metric_uu) << '\n';
    std::cerr << "metric_uu=\n" << metric_uu_func(position_u);
    std::cerr << "F_ll=\n" << F_ll;
    std::cerr << "F_ul=\n" << F_ul;
    std::cerr << "F_uu=\n" << F_uu;
    std::cerr << "F_ul*F_uu=\n" << F_ul * F_uu;
    std::cerr << "total_metric_uu=\n" << total_metric_uu << '\n';
    std::cerr << "total_metric_ll=\n" << total_metric_ll << '\n';
    std::cerr << "chr_lll=\n" << chr_lll;
  }
  Vector<T> result_l = Vector<T>();
  for (int k = 0; k < 4; ++k)
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        result_l[k] += vec_u[i] * chr_lll[k][i][j] * vec_u[j];
  for (int i = 0; i < 4; ++i)
    result_l[i] *= .5;  // Factor 1/2 from the Christoffel symbols.

  Vector<T> result_u = Vector<T>();
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      result_u[i] += total_metric_uu[i][j].value() * result_l[j];
  if (debug >= 3) {
    std::cerr << "result_u=" << result_u << '\n';
  }
  return result_u;
}

}  // namespace bhr

#endif
