#ifndef SPACETIME_H
#define SPACETIME_H

#include <bhr/tensor.h>

#define DERIVED (static_cast<const Derived *>(this))

// https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
template <typename Derived>
class SpacetimeBase {
 public:
  template <typename Vector>
  inline auto dot(const Vector &position,
                  const Vector &A,
                  const Vector &B) const {
    typedef typename Vector::value_type T;
    Matrix4<T> metric = DERIVED->get_metric_ll(position);

    T result = T();
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        result += A[i] * metric[i][j] * B[j];
    return result;
  }

  template <typename T, template <typename> class Vector>
  inline Vector<T> raise_index(const Vector<T> &position,
                                 const Vector<T> &vec) const {
    Matrix4<T> metric = DERIVED->get_metric_uu(position);

    Vector<T> result = Vector<T>();
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        result[i] += metric[i][j] * vec[j];
    return result;
  }

  template <typename T, template <typename> class Vector>
  inline Christoffel<T> get_christoffel_lll(
      const Vector<T> &position) const {
    typedef first_partial_derivatives<T, 4> fpds;
    Vector<fpds> ad_position{
      fpds(position[0], (T)1, (T)0, (T)0, (T)0),
      fpds(position[1], (T)0, (T)1, (T)0, (T)0),
      fpds(position[2], (T)0, (T)0, (T)1, (T)0),
      fpds(position[3], (T)0, (T)0, (T)0, (T)1)
    };

    Matrix4<fpds> ad_metric = DERIVED->get_metric_ll(ad_position);
    Christoffel<T> christoffel_lll;
    for (int k = 0; k < 4; ++k)
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
          christoffel_lll[k][i][j] = (ad_metric[k][j].first(i)
                                    + ad_metric[k][i].first(j)
                                    - ad_metric[i][j].first(k)) / 2;
        }
    return christoffel_lll;
  }

  template <typename T, template <typename> class Vector>
  inline Christoffel<T> get_christoffel_ull(
      const Vector<T> &position) const {
    Christoffel<T> christoffel_lll = DERIVED->get_christoffel_lll(position);
    Christoffel<T> christoffel_ull;
    Matrix4<T> metric_uu = DERIVED->get_metric_uu(position);
    for (int k = 0; k < 4; ++k)
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
          T tmp = T();
          for (int l = 0; l < 4; ++l)
            tmp += metric_uu[k][l] * christoffel_lll[l][i][j];
          christoffel_ull[k][i][j] = tmp;
        }
    return christoffel_ull;
  }

  template <typename T, template <typename> class Vector>
  inline Vector<T> geodesic_acceleration(
      const Vector<T> &position,
      const Vector<T> &direction) const {
    Christoffel<T> christoffel_lll = DERIVED->get_christoffel_lll(position);
    Vector<T> result_l = Vector<T>();
    for (int k = 0; k < 4; ++k)
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
          result_l[k] += direction[i] * christoffel_lll[k][i][j] * direction[j];

    return DERIVED->raise_index(position, result_l);
  }

  template <typename T, template <typename> class Vector>
  inline std::pair<Vector<T>, Vector<T>>
  geodesic_acceleration__parallel_transport(
      const Vector<T> &position,
      const Vector<T> &direction,
      const Vector<T> &additional) const {
    Christoffel<T> christoffel_lll = DERIVED->get_christoffel_lll(position);
    Vector<T> result_l = Vector<T>();
    for (int k = 0; k < 4; ++k)
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
          result_l[k] += direction[i] * christoffel_lll[k][i][j] * direction[j];

    Vector<T> result2_l = Vector<T>();
    for (int k = 0; k < 4; ++k)
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
          result2_l[k] +=
              direction[i] * christoffel_lll[k][i][j] * additional[j];
        }

    return std::make_pair(
        DERIVED->raise_index(position, result_l),
        DERIVED->raise_index(position, result2_l)
    );
  }
};


// https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
template <typename Derived>
class BlackHoleBase : public SpacetimeBase<Derived> {
 public:
  template <class Vector, typename T = typename Vector::value_type>
  inline bool is_in_black_hole(const Vector &pos,
                               const T &epsilon = T(1e-5),
                               const T &factor = T(1)) const {
    return pos.get_r() < factor * DERIVED->black_hole_radius() + epsilon;
  }

  template <typename Vector>
  inline bool is_in_near_flat(const Vector &position,
                              real_t flat_measure) const {
    return position.get_r() > DERIVED->black_hole_radius() / flat_measure;
  }
};

template <typename T, typename Spacetime, template <typename> class _Vector4>
inline T guess_null_geodesic_t_coord(const Spacetime &spacetime,
                                      const _Vector4<T> &position,
                                      const _Vector4<T> &direction) {
  /* Determine direction.t. */
  Matrix4<T> metric = spacetime.get_metric_ll(position);
  T a = metric[0][0];
  T b = metric[0][1] * direction[1]
       + metric[0][2] * direction[2]
       + metric[0][3] * direction[3];
  T c = T();

  // Solve a t^2 + 2 b t + c == 0.
  for (int i = 1; i < 4; ++i)
    for (int j = 1; j < 4; ++j)
      c += metric[i][j] * direction[i] * direction[j];

  T D = sqr(b) - a * c;
  if (D < 0) {
    fprintf(stderr, "D not >= 0!!! D = %lf\n", D);
    exit(1);
  }

  /* Assuming g_tt < 0 and direction.t > 0. */
  return (-std::sqrt(D) - b) / a;
}

#undef DERIVED
#endif
