#ifndef RAYTRACER_H
#define RAYTRACER_H

#include <bhr/coordinate.h>
#include <bhr/integration.h>
#include <bhr/physical_constants.h>
#include <bhr/utility.h>

namespace bhr {

class Camera {
 public:
  bool horizontal_flip = false;

  virtual ~Camera() {}
  virtual std::pair<Vector3, Vector3> ray(real_t x, real_t y) const = 0;
};

class OrthographicProjectionCamera : public Camera {
 private:
  Vector3 view, right;
 public:
  Vector3 eye, center, up;
  real_t aspect;
  real_t ver_range;

  OrthographicProjectionCamera() {}
  OrthographicProjectionCamera(const Vector3 &_eye,
                               const Vector3 &_center,
                               const Vector3 &_up,
                               real_t _aspect)
        : eye(_eye), center(_center), up(_up), aspect(_aspect) {
    _recalculate();
  }
  virtual ~OrthographicProjectionCamera() {}

  void _recalculate(void) {
    up.normalize();

    view = center - eye;
    view.normalize();

    right = cross(view, up);
    right.normalize();

    up = cross(right, view);
  }

  void set_aspect(real_t _aspect) {
    aspect = _aspect;
  }

  virtual std::pair<Vector3, Vector3> ray(real_t x, real_t y) const {
    if (horizontal_flip)
      x = -x;
    Vector3 result(eye);
    result += (.5 * ver_range * aspect * x) * right;
    result -= (.5 * ver_range * y) * up;
    return std::make_pair(result, view);
  }
};

class ProjectionCamera : public Camera {
 private:
  Vector3 view, right;
  real_t fovy, fovy_tan;
 public:
  Vector3 eye, center, up;
  real_t aspect;

  ProjectionCamera() {}
  ProjectionCamera(const Vector3 &_eye, const Vector3 &_center,
      const Vector3 &_up, real_t _aspect)
          : eye(_eye), center(_center), up(_up), aspect(_aspect) {
    _recalculate();
    set_fov(60);
  }
  virtual ~ProjectionCamera() {}

  void _recalculate(void) {
    up.normalize();

    view = center - eye;
    view.normalize();

    right = cross(view, up);
    right.normalize();

    up = cross(right, view);
  }

  void set_fov(real_t _fovy) {
    fovy = _fovy / 180 * M_PI / 2.0;
    fovy_tan = tan(fovy);
  }

  void set_aspect(real_t _aspect) {
    aspect = _aspect;
  }

  virtual std::pair<Vector3, Vector3> ray(real_t x, real_t y) const {
    if (horizontal_flip)
      x = -x;
    Vector3 result(view);
    result += (fovy_tan * aspect * x) * right;
    result -= (fovy_tan * y) * up;
    return std::make_pair(eye, result.normalize());
  }
};



template <typename BasicGeodesicState, typename Spacetime, typename Field>
std::pair<real_t, real_t> advance_geodesic__RKF45(
    const Spacetime &spacetime,
    const Field &field,
    const real_t min_h,
    real_t h,
    const real_t max_h,
    const real_t epsilon,
    const BasicGeodesicState &basic,
    BasicGeodesicState * const output) {
  /* basic must not be the same object as *output. */

  constexpr real_t SAFETY = real_t(0.84);
  constexpr real_t SAFETY_INC = real_t(1.2);
  constexpr real_t SAFETY_INC_MAX = real_t(3.0);

  auto RHS = [&spacetime, &field](const BasicGeodesicState &st) {
    return st.integration_step(spacetime, field);
  };

  int limit = 0;
  for (;;) {
    if (debug >= 3) {
      fprintf(stderr, "  h=%lg [%lg %lg]\n",
          (double)h, (double)min_h, (double)max_h);
    }
    const real_t R = integration_step__RKF45(RHS, h, basic, output);
    const real_t delta = SAFETY * std::pow(epsilon / R, 0.25);
    if (debug >= 3) {
      fprintf(stderr, "  h=%lg [%lg %lg]  R=%lg epsilon=%lg delta=%lg  ",
          (double)h, (double)min_h, (double)max_h,
          (double)R, (double)epsilon, (double)delta);
      std::cerr << output->position << "  " << output->direction << "\n";
    }
    if (R <= epsilon) {
      double new_h = delta > SAFETY_INC_MAX
          ? std::min(h * SAFETY_INC_MAX, max_h)
          : (delta > SAFETY_INC ? std::min(h * delta / SAFETY_INC, max_h) : h);
      return std::make_pair(h, new_h);
    } else if (h == min_h) {
      return std::make_pair(h, h);
    }
    h = clamp(h * delta, min_h, max_h);
    if (++limit == 20) {
      fprintf(stderr, "  h=%lg [%lg %lg]  R=%lg epsilon=%lg delta=%lg  ",
          (double)h, (double)min_h, (double)max_h,
          (double)R, (double)epsilon, (double)delta);
      std::cerr << output->position << "  " << output->direction << "\n";
      fprintf(stderr, "advance_geodesic__RKF45 step count limit exceeded.\n");
      exit(1);
    }
  }
}




template <typename FullGeodesicData, typename Coord,
          typename Spacetime, typename Field, typename Func>
auto generate_geodesic(
    const Spacetime &spacetime,
    const Field &field,
    Func break_condition_func,
    Coord position,
    Coord direction,
    int N,
    real_t dlambda,
    FullGeodesicData *output)
        -> decltype(break_condition_func(position, direction)) {

  typedef decltype(break_condition_func(position, direction)) return_type;
  // typedef typename Coord::value_type T;

  const real_t dlambda_min = dlambda / 100000;
  const real_t dlambda_max = dlambda * 100;

  typename FullGeodesicData::basic_type yn(position, direction), yn1;

  for (int n = 0; n < N; ++n) {
    if (yn.position.get_r() < NEUTRON_STAR_r * 10) {
      // dlambda = std::max(dlambda / 3, dlambda_min);
      double x = (yn.position.get_r() / NEUTRON_STAR_r - 1) / 9;
      // x = sqr(x);
      // double new_dlambda = dlambda_min * std::pow(
      //     dlambda_max / dlambda_min, x);
      double new_dlambda = dlambda_min * (1 - x) + dlambda_max * x;
      dlambda = std::min(dlambda, new_dlambda);
    }

    if (yn.position.is_spherical) {
      double theta = std::abs(Mod(yn.position.spherical_part(
              spacetime.coord_system_parameters(yn.position)).theta, M_PI));
      theta = std::min(theta, M_PI - theta);
      if (theta < 0.10) {
        double new_dlambda = dlambda_min * std::pow(
            dlambda_max / dlambda_min, .5 + 100 * theta);
        dlambda = std::min(dlambda, new_dlambda);
      }
    }

    std::pair<real_t, real_t> _dlambda = advance_geodesic__RKF45(
        spacetime,
        field,
        dlambda_min,
        dlambda,
        dlambda_max,
        1e-7,
        yn,
        &yn1);

    if (debug >= 10) {
      fprintf(stderr, "%7d ", n);
      std::cerr << dlambda << "  " << yn.position << "  " << yn.direction << '\n';
    }
    auto break_result = break_condition_func(yn1.position, yn1.direction);
    if (break_result) {
      if (debug >= 10) fprintf(stderr, "%d\n", n);
      output->basic = break_result != DEAD_BLACK_HOLE ? yn1 : yn;
      output->extra.finish(break_result, n);
      if (N >= 3) assert(n >= 3);
      return break_result;
    }

#if RENDER_DISK
    {
      const auto z0 = yn.position.get_z();
      const auto z1 = yn1.position.get_z();
      const auto r0 = yn.position.get_r();
      const auto r1 = yn1.position.get_r();
#if RENDER_DISK == DISK_SHAKURA
      const double h = shakura_sunyaev_height(r0);  // Hack.
#else
      constexpr double h = 0.0;
#endif
      int region0 = z0 >= h ? 1 : (z0 <= -h ? 3 : 2);
      int region1 = z1 >= h ? 1 : (z1 <= -h ? 3 : 2);
      if (region0 == 2 && INNER_RADIUS < r0 && r0 < OUTER_RADIUS) region0 |= 4;
      if (region1 == 2 && INNER_RADIUS < r1 && r1 < OUTER_RADIUS) region1 |= 4;
      if ((region0 != region1 || region0 == 6 || region1 == 6)
          && 0.5 * INNER_RADIUS < r0 && r0 < 1.5 * OUTER_RADIUS) {
        if (dlambda == dlambda_min) {
          if (INNER_RADIUS < r0 && r0 < OUTER_RADIUS) {
            if (N % 2) fprintf(stderr, "%d\n", n);
            output->basic = yn;
            output->extra.finish(DEAD_DISK, n);
            return DEAD_DISK;
          }
        } else if (0.7 * INNER_RADIUS < r0 && r0 < 1.3 * OUTER_RADIUS) {
          dlambda = std::max(dlambda_min, dlambda / 2);
          continue;
        }
      }
    }
#endif

    dlambda = _dlambda.second;
    yn = yn1;
    output->extra.post_step(spacetime, yn, _dlambda.first);
  }

  if (debug >= 10) fprintf(stderr, "%d\n", N);
  output->basic = yn;
  output->extra.finish(return_type(), N);
  return return_type();
}

}  // namespace bhr

#endif
