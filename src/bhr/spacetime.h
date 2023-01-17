#ifndef SPACETIME_H
#define SPACETIME_H

#include <bhr/tensor.h>

#define DERIVED (static_cast<const _Derived *>(this))

// https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
template <typename _Derived>
class SpacetimeBase {
 public:
  template <typename _Vector>
  inline auto dot(const _Vector &position,
                  const _Vector &A,
                  const _Vector &B) const {
    typedef typename _Vector::value_type _T;
    Matrix4<_T> metric = DERIVED->get_metric_ll(position);

    _T result = _T();
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        result += A[i] * metric[i][j] * B[j];
    return result;
  }

  template <typename _T, template <typename> class _Vector>
  inline _Vector<_T> raise_index(const _Vector<_T> &position,
                                 const _Vector<_T> &vec) const {
    Matrix4<_T> metric = DERIVED->get_metric_uu(position);

    _Vector<_T> result = _Vector<_T>();
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        result[i] += metric[i][j] * vec[j];
    return result;
  }

  template <typename _T, template <typename> class _Vector>
  inline Christoffel<_T> get_christoffel_lll(
      const _Vector<_T> &position) const {
    typedef first_partial_derivatives<_T, 4> fpds;
    _Vector<fpds> ad_position{
      fpds(position[0], (_T)1, (_T)0, (_T)0, (_T)0),
      fpds(position[1], (_T)0, (_T)1, (_T)0, (_T)0),
      fpds(position[2], (_T)0, (_T)0, (_T)1, (_T)0),
      fpds(position[3], (_T)0, (_T)0, (_T)0, (_T)1)
    };

    Matrix4<fpds> ad_metric = DERIVED->get_metric_ll(ad_position);
    Christoffel<_T> christoffel_lll;
    for (int k = 0; k < 4; ++k)
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
          christoffel_lll[k][i][j] = (ad_metric[k][j].first(i)
                                    + ad_metric[k][i].first(j)
                                    - ad_metric[i][j].first(k)) / 2;
        }
    return christoffel_lll;
  }

  template <typename _T, template <typename> class _Vector>
  inline Christoffel<_T> get_christoffel_ull(
      const _Vector<_T> &position) const {
    Christoffel<_T> christoffel_lll = DERIVED->get_christoffel_lll(position);
    Christoffel<_T> christoffel_ull;
    Matrix4<_T> metric_uu = DERIVED->get_metric_uu(position);
    for (int k = 0; k < 4; ++k)
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
          _T tmp = _T();
          for (int l = 0; l < 4; ++l)
            tmp += metric_uu[k][l] * christoffel_lll[l][i][j];
          christoffel_ull[k][i][j] = tmp;
        }
    return christoffel_ull;
  }

  template <typename _T, template <typename> class _Vector>
  inline _Vector<_T> geodesic_acceleration(
      const _Vector<_T> &position,
      const _Vector<_T> &direction) const {
    Christoffel<_T> christoffel_lll = DERIVED->get_christoffel_lll(position);
    _Vector<_T> result_l = _Vector<_T>();
    for (int k = 0; k < 4; ++k)
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
          result_l[k] += direction[i] * christoffel_lll[k][i][j] * direction[j];

    return DERIVED->raise_index(position, result_l);
  }

  template <typename _T, template <typename> class _Vector>
  inline std::pair<_Vector<_T>, _Vector<_T>>
  geodesic_acceleration__parallel_transport(
      const _Vector<_T> &position,
      const _Vector<_T> &direction,
      const _Vector<_T> &additional) const {
    Christoffel<_T> christoffel_lll = DERIVED->get_christoffel_lll(position);
    _Vector<_T> result_l = _Vector<_T>();
    for (int k = 0; k < 4; ++k)
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
          result_l[k] += direction[i] * christoffel_lll[k][i][j] * direction[j];

    _Vector<_T> result2_l = _Vector<_T>();
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
template <typename _Derived>
class BlackHoleBase : public SpacetimeBase<_Derived> {
 public:
  template <class _Vector, typename _T = typename _Vector::value_type>
  inline bool is_in_black_hole(const _Vector &pos,
                               const _T &epsilon = _T(1e-5),
                               const _T &factor = _T(1)) const {
    return pos.get_r() < factor * DERIVED->black_hole_radius() + epsilon;
  }

  template <typename _Vector>
  inline bool is_in_near_flat(const _Vector &position,
                              real_t flat_measure) const {
    return position.get_r() > DERIVED->black_hole_radius() / flat_measure;
  }
};

template <typename _T, typename _Spacetime, template <typename> class _Vector4>
inline _T guess_null_geodesic_t_coord(const _Spacetime &spacetime,
                                      const _Vector4<_T> &position,
                                      const _Vector4<_T> &direction) {
  /* Determine direction.t. */
  Matrix4<_T> metric = spacetime.get_metric_ll(position);
  _T a = metric[0][0];
  _T b = metric[0][1] * direction[1]
       + metric[0][2] * direction[2]
       + metric[0][3] * direction[3];
  _T c = _T();

  // Solve a t^2 + 2 b t + c == 0.
  for (int i = 1; i < 4; ++i)
    for (int j = 1; j < 4; ++j)
      c += metric[i][j] * direction[i] * direction[j];

  _T D = sqr(b) - a * c;
  if (D < 0) {
    fprintf(stderr, "D not >= 0!!! D = %lf\n", D);
    exit(1);
  }

  /* Assuming g_tt < 0 and direction.t > 0. */
  return (-std::sqrt(D) - b) / a;
}

#undef DERIVED
#endif
