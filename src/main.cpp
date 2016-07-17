#include <cassert>
#include <cmath>
#include <chrono>
#include <cstring>
#include <thread>

#include "base.h"
#include "utility.h"
#include "colorize.h"
#include "coordinate.h"

#include "field.h"
#include "flat.h"
#include "kerr.h"
// #include "reissner_nordstrom.h"
#include "schwarzschild.h"
#include "settings.h"
#include "tga.h"
#include "spectrum.h"
// #include "3rd/specrend.h"
#include "floatimage.h"

#include "../tests/tests.h"

int debug = 1;

// Hope this isn't reserved by something important.
#define GIVE_UP(...) { fprintf(stderr, __VA_ARGS__); return false; }

template <typename _FullGeodesicData,
          typename _Spacetime>
bool _generate_and_save_image(
    const _Spacetime &spacetime,
    const char *filename_float,
    const char *filename_image,
    Snapshot<_FullGeodesicData> *snapshot,
    RGBd *rgb,
#if SKY_ENABLED
    const Image &sky,
#endif
    const auto &disk_tex,
    int thread_count) {
#if RENDER_DISK && DISK_RELIEF_TEXTURE
  if (!load_disk_relief_texture())
    GIVE_UP("Disk texture loading failed. Aborting.\n");
#endif
  const int width = snapshot->width;
  const int height = snapshot->height;
  std::unique_ptr<RGBd[]> buffer_ptr;
  if (rgb == nullptr) {
    buffer_ptr.reset(rgb = new RGBd[width * height]);
#if SKY_ENABLED
    colorize(snapshot, spacetime, sky, disk_tex, rgb);
#else
    colorize(snapshot, spacetime, disk_tex, rgb);
#endif
  }

  if (filename_float && filename_float[0] != 0) {
    FILE *f = fopen(filename_float, "wb");
    if (f == nullptr)
      GIVE_UP("Error opening file %s!\n", filename_float);
    if (!save_floating_point_image(f, width, height, rgb)) {
      fclose(f);
      GIVE_UP("Error writing floating point image.\n");
    }
    fclose(f);
  }

  if (filename_image && filename_image[0] != 0) {
    std::unique_ptr<RGBA[]> rgba(new RGBA[width * height]);
    apply_color_filters(rgb, rgba.get(), width, height, thread_count);
    FILE *f = fopen(filename_image, "wb");
    if (f == nullptr)
      GIVE_UP("Error opening the file %s!\n", filename_image);
    fprintf(stderr, "Saving TGA to %s.\n", filename_image);
    if (!save_TGA(rgba.get(), width, height, f)) {
      fclose(f);
      GIVE_UP("Error saving TGA!\n");
    }
    fclose(f);

    // std::string filename_txt = std::string(filename_image) + ".csv";
    // if (!save_CSV(rgb, width, height, filename_txt.c_str()))
    //   GIVE_UP("Error saving CSV! Aborting!\n");
  }

  return true;
}

template <typename _FullGeodesicData>
bool _save_geodesics(
    const char *filename, Snapshot<_FullGeodesicData> *snapshot) {
  FILE *f = fopen(filename, "wb");
  if (f == nullptr) {
    fprintf(stderr, "Error opening file %s for writing!\n", filename);
    return false;
  }

  if (!snapshot->save(f)) {
    fprintf(stderr, "Error saving geodesics!\n");
    fclose(f);
    return false;
  }
  fclose(f);

  fprintf(stderr, "Saved geodesics to %s.\n", filename);
  return true;
}

template <typename _FullGeodesicData,
          typename _SpaceTime,
          typename _Field = Null>
bool generate_main(const Settings &S,
                   const _SpaceTime &spacetime,
                   const _Field &field = _Field()) {
  int width = S.width;
  int height = S.height;
  Camera *raytracer_camera;
  if (S.ortho) {
    OrthographicProjectionCamera *camera = new OrthographicProjectionCamera(
        S.camera_position_cart,
        S.camera_to_cart,
        S.camera_up_cart,
        real_t(width) / height);
    camera->ver_range = S.ver_range_km * UNIT_km;
    raytracer_camera = camera;
  } else {
    ProjectionCamera *camera = new ProjectionCamera(
        S.camera_position_cart,
        S.camera_to_cart,
        S.camera_up_cart,
        real_t(width) / height);
    camera->set_fov(S.fovy);
    raytracer_camera = camera;
  }
  raytracer_camera->horizontal_flip = S.horizontal_flip;
  if (debug) std::cerr << "Camera position=" << S.camera_position_cart << '\n';

#if RENDER_DISK && DISK_POLARIZATION
  // assert(!S.ortho);
  if (S.ortho) {
    // Using camera->eye instead of the infinitely distant observer.
    fprintf(stderr, "WARNING: Ortho not fully supported!\n");
  }
  typedef typename _FullGeodesicData::basic_type::coord_type _Coord;
  typedef typename _Coord::value_type _T;
  Vector3 eye =
      S.ortho ? ((OrthographicProjectionCamera *)raytracer_camera)->eye
              : ((ProjectionCamera *)raytracer_camera)->eye;
  CartesianVector4<_T> camera_position_cart(extend_vector((_T)0, eye));
  _Coord camera_position;
  convert_point(Null(), camera_position_cart,
                spacetime.coord_system_parameters(_Coord()), camera_position);
# if RENDER_DISK == DISK_KERTAP
    KERTAPDiskTexture<_Coord> disk_tex(camera_position);
# else
# error What kind of a disk with polarization enabled?
# endif
#elif RENDER_DISK == DISK_DUMMY
  DummyDiskTexture disk_tex;
#elif RENDER_DISK == DISK_SHAKURA
    ShakuraSunyaevDisk disk_tex(
        BLACK_HOLE_M,
        1e-14 * PHY_M_Sun / UNIT_y,
        _KERR_ISCO,
        1e15 * UNIT_Hz  //  4eV
        );
#elif RENDER_DISK
#error What kind of a disk?
#else
  DummyTexture disk_tex;
#endif

  real_t dlambda = S.dlambda;
  auto dlambda_func = [dlambda](const auto &/* pos */,
                                const auto &/* dir */) {
    return dlambda;
    // real_t dot = pos[1] * dir[1] + pos[2] * dir[2];
    // real_t rr = sqr(pos[1]) + sqr(pos[2]);
    // real_t pp = sqr(dir[1]) + sqr(dir[2]);
    // real_t dd = rr - sqr(dot) / pp;
    // return dlambda * (1 - 0.9 * exp(-dd / 0.001));
  };

  std::unique_ptr<Snapshot<_FullGeodesicData>> snapshot(S.recursive
      ? (Snapshot<_FullGeodesicData> *)new SnapshotRecursive<_FullGeodesicData>(width, height)
      : (Snapshot<_FullGeodesicData> *)new SnapshotMatrix<_FullGeodesicData>(width, height));
  RGBd *float_image = nullptr;
  std::string auto_geodesics_filename = S.get_auto_geodesics_filename();
  std::string output_float = S.get_float_filename();
  std::string output_image = S.get_image_filename();
  bool should_render = true;

  if (!S.cache && !output_float.empty() && !output_image.empty())
    GIVE_UP("No output specified!\n");

  if (S.input_filename.empty()) {
    if (S.cache) {
      FILE *f = fopen(auto_geodesics_filename.c_str(), "rb");
      if (f == nullptr) {
        fprintf(stderr, "Geodesics cache file %s not found, rendering.\n",
            auto_geodesics_filename.c_str());
      } else {
        if (!snapshot->load(f)) {
          fclose(f);
          GIVE_UP("Snapshot loading failed! Aborting!\n");
        }
        fprintf(stderr, "Successfully read geodesics cache file %s.\n",
            auto_geodesics_filename.c_str());
        fclose(f);
        should_render = false;
      }
    } else {
      fprintf(stderr, "Cache disabled.\n");
    }
  } else if (ends_with(S.input_filename, ".flt")) {
    FILE *f = fopen(S.input_filename.c_str(), "rb");
    if (f == nullptr) {
      GIVE_UP("Floating point image file %s not found! Aborting!\n",
          S.input_filename.c_str());
    }
    float_image = load_floating_point_image(f, width, height);
    fclose(f);
    if (float_image == nullptr)
      GIVE_UP("Error reading floating point image! Aborting!\n");
    fprintf(stderr, "Successfully read floating point image %s.\n",
        S.input_filename.c_str());
    should_render = false;
  } else if (ends_with(S.input_filename, ".geo")) {
    FILE *f = fopen(S.input_filename.c_str(), "rb");
    if (f == nullptr)
      GIVE_UP("Geodesics input file %s not found! Aborting!\n",
          S.input_filename.c_str());
    if (!snapshot->load(f)) {
      fclose(f);
      GIVE_UP("Loading failed, aborting!\n");
    }
    fprintf(stderr, "Successfully read geodesics input file %s.\n",
        S.input_filename.c_str());
    should_render = false;
  } else {
    GIVE_UP("Unrecognized input file extension %s. Aborting!\n",
        S.input_filename.c_str());
  }

  if (should_render) {
#define RECURSIVE 0
    if (S.recursive) {
#if RECURSIVE
      generate_image_recursive(spacetime, field, raytracer_camera, dlambda_func,
          width, height,
          reinterpret_cast<SnapshotRecursive<_FullGeodesicData> *>(snapshot.get()),
          S.threads, S.max_extra_recursive_depth);
#else
      assert(0 == 1); // Disabled to speed up compilation.
#endif
    } else {
#if RECURSIVE
      assert(0 == 1); // Disabled to speed up compilation.
#else
      generate_image(
          spacetime, field, raytracer_camera, dlambda_func,
          width, height,
          reinterpret_cast<SnapshotMatrix<_FullGeodesicData> *>(snapshot.get()),
          S.threads);
#endif
#undef RECURSIVE
    }
    if (S.cache)
      if (!_save_geodesics(auto_geodesics_filename.c_str(), snapshot.get()))
        GIVE_UP("Error saving geodesics, aborting!\n");
    snapshot->save_extra();
  }

  if (!output_float.empty() || !output_image.empty()) {
#if SKY_ENABLED
    const Image *sky = load_image_and_cache(S.sky_image.c_str());
    if (sky == nullptr)
      GIVE_UP("Error loading sky texture, aborting!\n", S.sky_image.c_str());
#endif

    if (!_generate_and_save_image(
        spacetime,
        output_float.c_str(),
        output_image.c_str(),
        snapshot.get(),
        float_image,
#if SKY_ENABLED
        *sky,
#endif
        disk_tex,
        S.threads)) {
      GIVE_UP("Error generating or saving the image! Aborting!\n");
    }
  }

  return true;
}

// template <typename _Coord> using __FullGeodesicData =
//     FullGeodesicData<
//       BasicGeodesicState<_Coord>,
//       GeodesicExtraBase<
//         GeodesicExtra__ParallelTransport<BasicGeodesicState<_Coord>>,
//         GeodesicExtra__MinR,
//         GeodesicExtra__Steps,
//         GeodesicExtra__DeadReason
//       >
//     >;

template <typename _Coord>
struct _BasicGeodesicState
    : BasicGeodesicState<_Coord, Geodesic__ParallelTransport<_Coord>> {
  typedef BasicGeodesicState<_Coord, Geodesic__ParallelTransport<_Coord>> _Base;

  using _Base::_Base;

  friend inline auto numerical_sqr_distance(
      const _BasicGeodesicState &A, const _BasicGeodesicState &B) {
    return numerical_sqr_distance(_Base(A), _Base(B));
  }

  friend inline void mult(_BasicGeodesicState *self, const auto &c) {
    self->_Base::_mult(c);
  }
  friend inline void mult_add(
      _BasicGeodesicState *self, const auto &c, const _BasicGeodesicState &A) {
    self->_Base::_mult_add(c, A);
  }
  friend inline void set_and_mult_add(
      _BasicGeodesicState *self, const _BasicGeodesicState &A, const auto &c,
      const _BasicGeodesicState &B) {
    self->_Base::_set_and_mult_add(A, c, B);
  }

  _BasicGeodesicState integration_step(const auto &spacetime,
                                       const auto &/* field */) const {
    _BasicGeodesicState result;
    typedef typename _BasicGeodesicState::coord_type::value_type _T;
    auto christoffel_ull = spacetime.get_christoffel_ull(this->position);
    for (int k = 0; k < 4; ++k) {
      _T tmp = _T();
      for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
          tmp += this->direction[i]
              * christoffel_ull[k][i][j]
              * this->direction[j];
      result.direction[k] = tmp;
    }
    result.position = this->direction;
    mult(&result.position, -1);
    result.Geodesic__ParallelTransport<_Coord>::__integration_step__impl(
        christoffel_ull, *this);
    return result;
  }
};

#if MAGNETIC_FIELD
template <typename _Coord> using __FullGeodesicData =
    FullGeodesicData<
      BasicGeodesicState<_Coord>,
      GeodesicExtraBase<
        GeodesicExtra__MinR,
        GeodesicExtra__Steps,
        GeodesicExtra__Debug,
        GeodesicExtra__DeadReason
      >
    >;
#elif DISK_POLARIZATION
// TODO #if DISK_POLARIZATION
template <typename _Coord> using ___FullGeodesicData =
    FullGeodesicData<
      _BasicGeodesicState<_Coord>,
      GeodesicExtraBase<
        GeodesicExtra__MinR,
        GeodesicExtra__Steps,
//        GeodesicExtra__SavedLambdas,
//        GeodesicExtra__Debug,
        GeodesicExtra__DeadReason
      >
    >;

template <typename _Coord>
struct __FullGeodesicData : ___FullGeodesicData<_Coord> {
  using ___FullGeodesicData<_Coord>::___FullGeodesicData;
  typedef FullGeodesicData<
      BasicGeodesicState<_Coord>,
      GeodesicExtraBase<
        GeodesicExtra__MinR,
        GeodesicExtra__Steps,
        GeodesicExtra__DeadReason
      > > simple_version;

  __FullGeodesicData() {}
  __FullGeodesicData(const simple_version &simple) {
    this->basic.position = simple.basic.position;
    this->basic.direction = simple.basic.direction;
    this->extra = simple.extra;
  }
};
#else
template <typename _Coord> using __FullGeodesicData =
    FullGeodesicData<
      BasicGeodesicState<_Coord>,
      GeodesicExtraBase<
        GeodesicExtra__MinR,
        GeodesicExtra__Steps,
//        GeodesicExtra__Debug,
        GeodesicExtra__DeadReason
      >
    >;
#endif

int main(int argc, char **argv) {
#if CHECK_KERR
  check_kerr();
  return 0;
#endif

#if RENDER_DISK
  calculate_spectrum_colors();
#endif
  // for (int i = 0; i < 100; ++i) {
  //   double T = i * 200;
  //   RGBd rgb = get_black_body_color(T);
  //   bool changes = constrain_rgb(&rgb);

  //   printf("T=%lf  %lf %lf %lf  %s\n", T, rgb[0], rgb[1], rgb[2],
  //       changes ? "(*)" : "");
  // }
  // return 0;

  Settings S;
  S.set_default();
  if (!S.read_cli(argc, argv)) return 1;
  if (!S.check()) return 1;
  debug = S.debug;

  if (debug) debug_units();

#if PREPROCESS_LAMBDAS
  if (!QEDCache::lambda_preprocess(S.threads))
    return 1;
#if CHECK_LAMBDA_PRECISION
  QEDCache::check_lambda_interpolation_precision();
  return 0;
#endif
#if GENERATE_LAMBDAS
  QED::generate_lambdas();
  return 0;
#endif
#elif CHECK_LAMBDA_PRECISION
#error Cannot check lambda interpolation precision when cache is disabled.
#elif GENERATE_LAMBDAS
#error Cannot generate lambdas when cache is disabled.
#endif

#if TESTS_ENABLED
  if (!test_all())
    return 1;
  return 0;
#endif

  auto start_time = std::chrono::system_clock::now();
  bool ok = false;

//   switch (S.type) {
#if 0
    // case SPACETIME_FLAT:
      FlatSpacetime spacetime;
      // FlatDipole field(S.field[0] * UNIT_T,
      //                  S.field[1] * M_PI / 180,
      //                  S.field[2] * M_PI / 180);
      Null field;
      // typedef CartesianVector4<real_t> _Coord;
      typedef SphericalVector4<real_t> _Coord;
      // ok = generate_main<__FullGeodesicData<CartesianVector4<real_t>>>(
      //     S, FlatSpacetime(), field);
      // ok = generate_main<SphericalVector4<real_t>>(S, FlatSpacetime());
    // break;
#endif
#if PREDEFINED_PARAMS
#if 1
    // case SPACETIME_KERR:
    typedef BoyerLindquistVector4<real_t> _Coord;
    KerrSpacetime spacetime;
    Null field;
    if (debug) {
      fprintf(stderr, "(parse debug script scale) _BLACK_HOLE_r_S = %lg\n",
        sqrt(_BLACK_HOLE_r_S * spacetime.black_hole_radius()));
    }
    //  ok = generate_main<__FullGeodesicData<BoyerLindquistVector4<real_t>>>(
    //      S, KerrSpacetime());
    // break;
#endif
#if 0
//     case SPACETIME_SCHWARZSCHILD:
      SchwarzschildSpacetime spacetime;
      // FlatSpacetime spacetime;
      SchwarzschildDipole field(spacetime, S.field[0] * UNIT_T);
      // Null field;
      typedef SphericalVector4<real_t> _Coord;
      // _Coord tmp = _Coord{0., NEUTRON_STAR_r, M_PI / 2, 0.};
      // std::cerr << "F_ll [Tm otprilike]\n" << (1 / UNIT_T / UNIT_m) * field.get_F_ll(tmp);
      // std::cerr << "F=" << field._calc_F(tmp) / sqr(UNIT_T) << " T^2\n";
      // debug = 3;
      // geodesic_acceleration__magnetic_field(
      //   [&](auto position_u) { return spacetime.get_metric_ll(position_u); },
      //   [&](auto position_u) { return spacetime.get_metric_uu(position_u); },
      //   [&](auto position_u) { return field.get_F_ll(position_u); },
      //   [&](auto F, auto G) { return EH::lagrangian_real(F, G); },
      //   tmp,
      //   tmp  // direction
      // );
      // exit(1);
//       break;
#endif
    // case SPACETIME_REISSNER_NORDSTROM:
    //   ok = generate_main<SphericalVector4<real_t>>(
    //       S, ReissnerNordstromSpacetime());
    //   break;
    // case SPACETIME_SCHWARZSCHILD:
    //   ok = generate_main<SphericalVector4<real_t>>(
    //       S, SchwarzschildSpacetime(S.M));
    //   break;
    // case SPACETIME_KERR:
    //   ok = generate_main<BoyerLindquistVector4<real_t>>(
    //       S, KerrSpacetime(S.M, S.a));
    //   break;
    // case SPACETIME_REISSNER_NORDSTROM:
    //   ok = generate_main<SphericalVector4<real_t>>(
    //       S, ReissnerNordstromSpacetime(S.M, S.Q));
    //   break;
#endif
//    default:
//      fprintf(stderr, "Unrecognized spacetime type %d\n", S.type);
//      return 1;
//  };


  // for (double x = 1.01; x <= 10.; x += 0.1) {
  //   double r = spacetime.black_hole_radius() * x;
  //   double theta = M_PI / 2;
  //   double phi = 0;
  //   SphericalVector4<double> position{0, r, theta, phi};
  //   double F = field._calc_F(position);
  //   double B = std::sqrt(2 * F);
  //   fprintf(stderr, "x=%lf B=%lg T\n", x, B / UNIT_T);
  // }
  // fprintf(stderr, "NEUTRON_STAR_r / BH_r = %lg\n",
  //   NEUTRON_STAR_r / spacetime.black_hole_radius());
  // return 0;

  if (S.height <= 0) {

    // FlatSpacetime spacetime;
    // FlatDipole field(S.field[0] * UNIT_T,
    //                  S.field[1] * M_PI / 180,
    //                  S.field[2] * M_PI / 180);
    // SchwarzschildSpacetime spacetime;
    // Null field;
    OrthographicProjectionCamera *camera = new OrthographicProjectionCamera(
        S.camera_position_cart,
        S.camera_to_cart,
        S.camera_up_cart,
        1.);
    camera->ver_range = S.ver_range_km * UNIT_km;

    const int N = std::max(1, S.width);
    const double dlambda = S.dlambda;
    for (int i = (-S.height); i < N; ++i) {
      double x = (double)i / (N - 1);
      double y = 0;
      __FullGeodesicData<_Coord> result;
      integrate_single_geodesic(
          spacetime,
          field,
          camera,
          [dlambda](auto, auto){ return dlambda; },
          x, y,
          &result);
      auto &position = result.basic.position;
      auto &direction = result.basic.direction;
      // Vector3 pos3{position[1], position[2], position[3]};
      // Vector3 dir3{direction[1], direction[2], direction[3]};
      // std::cerr << position << '\n';
      // double cos_theta = dir3.dot(pos3) / (pos3.length() * dir3.length());
      // std::cerr << "cos_theta=" << cos_theta << '\n';
      // std::cerr << "direction.length=" << dir3.length() << '\n';
      // position += (1000 - position.get_r()) / cos_theta / dir3.length() * direction;

      if (S.height == 0) {
        std::cerr << position << '\t' << direction << '\n';
        auto spherical = result.basic.position.spherical_part(Null());
        std::cerr << spherical.theta << "   " << spherical.phi << '\n';
      }

      // fprintf(stderr, "%.9lf %.9lf\n",
      //     position[0],
      //     position.get_r());
    }
    return 0;
  } else {
    ok = generate_main<__FullGeodesicData<_Coord>>(S, spacetime, field);
  }

  std::chrono::duration<double> time_delta =
      std::chrono::system_clock::now() - start_time;
  fprintf(stderr, "Finished in %.2lfs.\n\n", time_delta.count());

  return ok ? 0 : 1;
}
