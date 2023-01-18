#ifndef GEODESIC_H
#define GEODESIC_H

#include <bhr/utility.h>
#include <bhr/tensor.h>
#include <bhr/matrix.h>

#if MAGNETIC_FIELD_FULL
#include <bhr/euler_heisenberg.h>
#endif

template <typename Spacetime, typename Field, typename Coord>
auto __direction(
    const Spacetime &spacetime,
    const Field &field,
    const Coord &position,
    const Coord &direction) {
#if MAGNETIC_FIELD_FULL
  return geodesic_acceleration__magnetic_field(
      [&](auto position_u) { return spacetime.get_metric_ll(position_u); },
      [&](auto position_u) { return spacetime.get_metric_uu(position_u); },
      [&](auto position_u) { return field.get_F_ll(position_u); },
      [&](auto F, auto G) { return EH::lagrangian_real(F, G); },
      position,
      direction
  );
#elif MAGNETIC_FIELD
  return geodesic_acceleration__magnetic_field__lowest_order(
      [&](auto position_u) { return spacetime.get_metric_ll(position_u); },
      [&](auto position_u) { return spacetime.get_metric_uu(position_u); },
      [&](auto position_u) { return field.get_potential_l(position_u); },
      position,
      direction
  );
#else
  (void)field;
  return spacetime.geodesic_acceleration(position, direction);
#endif
}

template <typename Spacetime, typename Coord>
auto __direction(
    const Spacetime &spacetime,
    const Null &/* field */,
    const Coord &position,
    const Coord &direction) {
  return spacetime.geodesic_acceleration(position, direction);
}

template <typename Coord, typename ...Extra>
struct BasicGeodesicState : Extra... {
  typedef Coord coord_type;
  Coord position;
  Coord direction;

  BasicGeodesicState() {}
  BasicGeodesicState(const Coord &_position, const Coord &_direction)
      : position(_position), direction(_direction) {}

  friend inline auto numerical_sqr_distance(
      const BasicGeodesicState &A, const BasicGeodesicState &B) {
    return numerical_sqr_distance(A.position, B.position)
           + (numerical_sqr_distance(A.direction, B.direction)
             + ... + numerical_sqr_distance(Extra(A), Extra(B)));
  }

  template <typename U>
  friend inline void mult(BasicGeodesicState *self, const U &c) {
    mult(&self->position, c);
    mult(&self->direction, c);
    (self->Extra::_mult(c), ...);
  }

  template <typename U>
  friend inline void mult_add(
      BasicGeodesicState *self, const U &c, const BasicGeodesicState &A) {
    mult_add(&self->position, c, A.position);
    mult_add(&self->direction, c, A.direction);
    (self->Extra::_mult_add(c, A), ...);
  }

  template <typename U>
  friend inline void set_and_mult_add(
      BasicGeodesicState *self,
      const BasicGeodesicState &A, const U &c, const BasicGeodesicState &B) {
    set_and_mult_add(&self->position, A.position, c, B.position);
    set_and_mult_add(&self->direction, A.direction, c, B.direction);
    (self->Extra::_set_and_mult_add(A, c, B), ...);
  }

  template <typename U>
  inline void _mult(const U &c) {
    mult(this, c);
  }

  template <typename U>
  inline void _mult_add(const U &c, const BasicGeodesicState &A) {
    mult_add(this, c, A);
  }

  template <typename U>
  inline void _set_and_mult_add(const BasicGeodesicState &A,
                                const U &c,
                                const BasicGeodesicState &B) {
    set_and_mult_add(this, A, c, B);
  }

  template <typename Spacetime, typename Field>
  inline BasicGeodesicState integration_step(
      const Spacetime &spacetime,
      const Field &field) const {
    BasicGeodesicState result;
    result.position = direction;
    result.direction = __direction(spacetime, field, position, direction);
    mult(&result.position, -1);
    (result.Extra::integration_step__impl(spacetime, *this), ...);
    return result;
  }
};


template <typename ...Extra>
struct GeodesicExtraBase : Extra... {
#if RENDER_DISK && RENDER_DISK != DISK_DUMMY
  // TEMPORARY: Used for disks.
  BoyerLindquistVector4<double> start_position;
  BoyerLindquistVector4<double> start_direction;
  // SphericalVector4<double> start_position;
  // SphericalVector4<double> start_direction;
#endif

  template <typename Spacetime, typename Basic>
  void post_step(const Spacetime &spacetime,
                 const Basic &basic,
                 const real_t &used_dlambda) {
    /* Call post_step__impl for each extra part. */
    (void)spacetime;
    (void)basic;
    (void)used_dlambda;
    (this->Extra::post_step__impl(spacetime, basic, used_dlambda), ...);
  }

  void finish(int dead_reason, int steps) {
    /* Call finish__impl for each extra part. */
    (void)dead_reason;
    (void)steps;
    (this->Extra::finish__impl(dead_reason, steps), ...);
  }

};

struct GeodesicExtra__Base {
  template <typename Spacetime, typename Basic>
  void post_step__impl(const Spacetime &/* spacetime */,
                       const Basic &/* basic */,
                       const real_t &/* used_dlambda */) {
    /* noop by default */
  }

  void finish__impl(int /* dead_reason */, int /* steps */) {
    /* noop by default */
  }
};

struct GeodesicExtra__MinR : GeodesicExtra__Base {
  real_t min_r = 1e+9;

  template <typename Spacetime, typename Basic>
  void post_step__impl(const Spacetime &/* spacetime */,
                       const Basic &basic,
                       const real_t &/* used_dlambda */) {
    min_r = std::min(min_r, basic.position.get_r());
  }
};

struct GeodesicExtra__Steps : GeodesicExtra__Base {
  int steps;

  void finish__impl(int /* dead_reason */, int steps) {
    this->steps = steps;
  }
};

struct GeodesicExtra__DeadReason : GeodesicExtra__Base {
  int dead_reason;

  void finish__impl(int dead_reason, int /* steps */) {
    this->dead_reason = dead_reason;
  }
};

struct GeodesicExtra__Debug : GeodesicExtra__Base {
  int n = 0;
  template <typename Spacetime, typename Basic>
  void post_step__impl(const Spacetime &spacetime,
                       const Basic &basic,
                       const real_t &used_dlambda) {
    // std::cerr << basic.position << '\n';
    // std::cerr << basic.parallel_transport_lu;
    if (debug <= 2) return;
    CartesianVector4<double> pos, dir;
    convert_point_and_diff(
        spacetime.coord_system_parameters(basic.position),
        basic.position,
        basic.direction,
        Null(),
        pos,
        dir);

    std::cerr << "n=" << n++
              << "\tused_dlambda=" << used_dlambda
              << "\t" << basic.position
              << "\t" << basic.direction
              << "\t<" << pos << ">"
              << "\t" << dir
              << "\n";
  }

  void finish__impl(int dead_reason, int steps) {
    if (debug <= 2) return;
    std::cerr << "dead_reason=" << dead_reason
              << "\tsteps=" << steps
              << "\n";
  }
};

struct GeodesicExtra__SavePath : GeodesicExtra__Base {
  std::vector<std::pair<CartesianVector4<double>,
                        CartesianVector4<double>>> path;
  template <typename Spacetime, typename Basic>
  void post_step__impl(const Spacetime &spacetime,
                       const Basic &basic,
                       const real_t &/* used_dlambda */) {
    CartesianVector4<double> pos, dir;
    convert_point_and_diff(
        spacetime.coord_system_parameters(basic.position),
        basic.position,
        basic.direction,
        Null(), pos, dir);
    path.emplace_back(pos, dir);
  }
};

struct GeodesicExtra__SavedLambdas : GeodesicExtra__Base {
  std::vector<double> dlambdas;
  template <typename Spacetime, typename Basic>
  void post_step__impl(const Spacetime &/* spacetime */,
                       const Basic &/* basic */,
                       const real_t &used_dlambda) {
    dlambdas.push_back(used_dlambda);
  }
};


template <typename Coord>
struct Geodesic__ParallelTransport {
  typedef typename Coord::value_type T;
  Matrix4<T> parallel_transport_lu;

  Geodesic__ParallelTransport() {
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < 4; ++j)
        parallel_transport_lu[i][j] = i == j ? 1 : 0;
  }

  template <typename Basic>
  void __integration_step__impl(const Christoffel<T> &christoffel_ull,
                                const Basic &basic) {
    for (int l = 0; l < 4; ++l) {
      for (int k = 0; k < 4; ++k) {
        T tmp = T();
        for (int i = 0; i < 4; ++i)
          for (int j = 0; j < 4; ++j) {
            tmp += basic.direction[i]
                 * christoffel_ull[k][i][j]
                 * basic.parallel_transport_lu[l][j];
          }
        parallel_transport_lu[l][k] = tmp;
      }
    }
  }

  template <typename Spacetime, typename Basic>
  void integration_step__impl(const Spacetime &spacetime, const Basic &basic) {
    Christoffel<T> christoffel_lll =
        spacetime.get_christoffel_lll(basic.position);
    for (int l = 0; l < 4; ++l) {
      Coord result;  // result_l
      for (int k = 0; k < 4; ++k) {
        T tmp = T();
        for (int i = 0; i < 4; ++i)
          for (int j = 0; j < 4; ++j) {
            tmp += basic.direction[i]
                 * christoffel_lll[k][i][j]
                 * basic.parallel_transport_lu[l][j];
          }
        result[k] = tmp;
      }

      // result_l -> result_u
      result = spacetime.raise_index(basic.position, result);
      // We do backwards raytracing, thus we have a + sign here.
      for (int k = 0; k < 4; ++k)
        parallel_transport_lu[l][k] = result[k];
    }
  }

  friend inline auto numerical_sqr_distance(
      const Geodesic__ParallelTransport &/* A */,
      const Geodesic__ParallelTransport &/* B */) {
    return 0.0;  // Exclude from the epsilon calculation.
    // return numerical_sqr_distance(A.parallel_transport_lu,
    //                               B.parallel_transport_lu);
  }

  template <typename U>
  inline void _mult(const U &c) {
    mult(&parallel_transport_lu, c);
  }

  template <typename U>
  inline void _mult_add(const U &c, const Geodesic__ParallelTransport &A) {
    mult_add(&parallel_transport_lu, c, A.parallel_transport_lu);
  }

  template <typename U>
  inline void _set_and_mult_add(const Geodesic__ParallelTransport &A,
                                const U &c,
                                const Geodesic__ParallelTransport &B) {
    set_and_mult_add(&parallel_transport_lu, A.parallel_transport_lu, c,
                     B.parallel_transport_lu);
  }
};


template <typename Basic, typename Extra>
struct FullGeodesicData {
  typedef Basic basic_type;
  typedef Extra extra_type;

  Basic basic;
  Extra extra;
};

template <typename Coord> using BasicFullGeodesicData =
    FullGeodesicData<BasicGeodesicState<Coord>, GeodesicExtraBase<>>;

#endif
