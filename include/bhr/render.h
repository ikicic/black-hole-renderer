#ifndef RENDER_H
#define RENDER_H

#include <bhr/disk.h>
#include <bhr/raytracer.h>
#include <bhr/geodesic.h>
#include <bhr/texture.h>
#include <bhr/spectrum.h>
#include <bhr/3rd/specrend.h>

#include <bhr/line.h>

#include <mutex>
#include <stack>
#include <thread>
#include <utility>

#define OBJECT_COLOR (RGBd{0., 1., 1.})
#if 1
#define SKY_COLOR (RGBd{0., 0., 0.})
#define ARROW_COLOR (RGBd{1., 1., 1.})
#define BLACK_HOLE_COLOR (RGBd{0., 0., 0.})
#else
#define SKY_COLOR (RGBd{1., 1., 1.})
#define ARROW_COLOR (RGBd{-1., -1., -1.})
// #define BLACK_HOLE_COLOR (RGBd{0.5, 0.5, 0.5})
#define BLACK_HOLE_COLOR (RGBd{0.7, 0.7, 0.7})
#endif

template <typename _FullGeodesicData>
class SnapshotMatrix : public Snapshot<_FullGeodesicData> {
 public:
  SnapshotMatrix(int width, int height)
      : Snapshot<_FullGeodesicData>(width, height), geodesics(nullptr) {}
  ~SnapshotMatrix() {
    delete []geodesics;
    geodesics = nullptr;
  }

  virtual bool load(FILE *f) {
    const int total = this->width * this->height;
    geodesics = new _FullGeodesicData[total];

    const int result = fread(geodesics, sizeof(geodesics[0]), total, f);
    if (result != total) {
      fprintf(stderr, "Error, read only %d/%d geodesics!\n", result, total);
      delete []geodesics;
      geodesics = nullptr;
      return false;
    }

    return true;
  }

  virtual bool save(FILE *f) const {
    const int total = this->width * this->height;
    if (fwrite(geodesics, sizeof(geodesics[0]) * total, 1, f) != 1) {
      fprintf(stderr, "Error writing geodesics!\n");
      return false;
    }
    return true;
  }

  virtual void __save_theta_phi(void) const {
    FILE *f2 = fopen("output/output_phi.csv", "w");
    FILE *f3 = fopen("output/output_theta.csv", "w");
    if (f2 == nullptr || f3 == nullptr) {
      fprintf(stderr, "__save_theta_phi can't open file.\n");
      exit(1);
    }
    for (int i = 0; i < this->height; ++i)
      for (int j = 0; j < this->width; ++j) {
        auto &geo = geodesics[i * this->width + j];
        auto spherical = geo.basic.position.spherical_part(Null());
        if (geo.extra.dead_reason == DEAD_BLACK_HOLE)
          spherical.phi = spherical.theta = 0;
#define FORMAT "%10.8lg"
        fprintf(f2, j == this->width - 1 ? FORMAT "\n" : FORMAT ",", spherical.phi);
        fprintf(f3, j == this->width - 1 ? FORMAT "\n" : FORMAT ",", spherical.theta);
#undef FORMAT
      }
    fclose(f2);
    fclose(f3);
  }

  _FullGeodesicData *geodesics;
};



void _generate_image_status_thread(std::vector<int> &finished_pixels_vec,
                                   int width,
                                   int height);

template <typename _Coord, typename _T, typename _Spacetime>
std::pair<_Coord, _Coord> cartesian3d_to_coord4d(
    const CartesianVector<_T, 3> &position_cart3,
    const CartesianVector<_T, 3> &direction_cart3,
    const _Spacetime &spacetime) {
  CartesianVector4<_T> position_cart4(extend_vector((_T)0, position_cart3));
  CartesianVector4<_T> direction_cart4(extend_vector((_T)0, direction_cart3));

  _Coord position, direction;
  convert_point_and_diff(
      Null(), position_cart4, direction_cart4,
      spacetime.coord_system_parameters(_Coord()), position, direction);
  direction[0] = guess_null_geodesic_t_coord(spacetime, position, direction);
  return std::make_pair(position, direction);
}


template <
    typename _Spacetime,
    typename _Field,
    typename _FullGeodesicData,
    typename DLambdaFunc>
int integrate_single_geodesic(
    const _Spacetime &spacetime,
    const _Field &field,
    const Camera *raytracer_camera,
    DLambdaFunc dlambda_func,
    real_t x,
    real_t y,
    _FullGeodesicData *output) {

  typedef typename _FullGeodesicData::basic_type basic_type;
  typedef typename basic_type::coord_type coord_type;
  // constexpr real_t flat_measure = 0.00001;

  // We do the backwards ray-tracing, so we're flipping the direction.
  coord_type position, direction;
  std::pair<Vector3, Vector3> ray = raytracer_camera->ray(x, y);
  std::tie(position, direction) = cartesian3d_to_coord4d<coord_type>(
      ray.first, -ray.second, spacetime);

  auto break_condition_func =
    [&spacetime](const coord_type &position, const coord_type &direction) {
      if (!std::isfinite(position) || !std::isfinite(direction))
        return DEAD_BLACK_HOLE;
      if (position.get_r() < NEUTRON_STAR_r)
        return DEAD_BLACK_HOLE;
      if (position.get_r() < _BLACK_HOLE_radius)
        return DEAD_BLACK_HOLE;
      if (spacetime.is_in_black_hole(position))
        return DEAD_BLACK_HOLE;
      if (position.get_r() > MAX_r) {
#if FAKE_SKY
        std::pair<double, double> grid = position.spherical_part(
              spacetime.coord_system_parameters(coord_type())).to_txty();
        int x = (int)(grid.first * 360);
        int y = (int)(grid.second * 180);

        return DEAD_SKY_TEX_OFFSET | (x << 10) | y;
#else
        return DEAD_FLAT;
#endif
      }
      return 0;
    };

  real_t dlambda = dlambda_func(position, direction);
  constexpr int N = 100000;
  if (N % 2) fprintf(stderr, "x=%.25lf y=%.25lf\n", x, y);
#if DISK_POLARIZATION
  typename _FullGeodesicData::simple_version simple;
  int result = generate_geodesic(
      spacetime, field, break_condition_func, position,
      direction, N, dlambda, &simple);
  if (result == DEAD_DISK) {
    result = generate_geodesic(spacetime, field, break_condition_func, position,
                               direction, N, dlambda, output);
  } else {
    *output = simple;
  }
#else
  int result = generate_geodesic(
      spacetime, field, break_condition_func,
      position, direction, N, dlambda, output);
#endif
  output->extra.start_position = position;
  output->extra.start_direction = direction;
  return result;
}

template <typename _Spacetime, typename _Field, typename _FullGeodesicData,
          typename _dlambdaFunc>
void _generate_image_worker_thread(
    std::stack<std::pair<int, int>> &tasks,
    std::mutex &tasks_mutex,
    int *finished_pixels,
    const _Spacetime &spacetime,
    const _Field &field,
    const Camera *raytracer_camera,
    const _dlambdaFunc &dlambda_func,
    int width,
    int height,
    _FullGeodesicData *output) {

  for (;;) {
    tasks_mutex.lock();
    if (tasks.empty()) {
      tasks_mutex.unlock();
      break;
    }
    std::pair<int, int> interval = tasks.top();
    tasks.pop();
    tasks_mutex.unlock();

    for (int i = interval.first; i <= interval.second; ++i) {
      for (int j = 0; j < width; ++j) {
        // i == 0, j == 0 --> x = -1, y = -1 --> Top left.
        double x = double(1 + 2 * j - width) / width;    // -1 <= x <= 1
        double y = double(1 + 2 * i - height) / height;  // -1 <= y <= 1
        integrate_single_geodesic(spacetime, field, raytracer_camera,
            dlambda_func, x, y, &output[i * width + j]);
        ++*finished_pixels;
      }
    }
  };
}

template <typename _FullGeodesicData, typename _Spacetime, typename _Field,
          typename _dlambdaFunc>
void generate_image(
    const _Spacetime &spacetime,
    const _Field &field,
    const Camera *raytracer_camera,
    const _dlambdaFunc &dlambda_func,
    int width,
    int height,
    SnapshotMatrix<_FullGeodesicData> *output_snapshot,
    int thread_count) {

  fprintf(stderr, "Rendering %dx%d image. Threads=%d\n",
      width, height, thread_count);

  delete []output_snapshot->geodesics;
  output_snapshot->geodesics = new _FullGeodesicData[width * height];

  std::stack<std::pair<int, int>> tasks;
  std::mutex tasks_mutex;

  constexpr int SPLIT = 1;
  for (int i = 0; i < (height + SPLIT - 1 )/ SPLIT; ++i)
    tasks.emplace(i * SPLIT, std::min(height - 1, (i + 1) * SPLIT - 1));

  std::vector<std::thread> threads;
  std::vector<int> finished_pixels_vec(thread_count);
  threads.emplace_back(
      _generate_image_status_thread,
      std::ref(finished_pixels_vec),
      width, height);
  for (int i = 0; i < thread_count; ++i) {
    threads.emplace_back(
        _generate_image_worker_thread<
            _Spacetime, _Field, _FullGeodesicData, _dlambdaFunc>,
        std::ref(tasks),
        std::ref(tasks_mutex),
        &finished_pixels_vec[i],
        std::ref(spacetime),
        std::ref(field),
        raytracer_camera,
        std::ref(dlambda_func),
        width,
        height,
        output_snapshot->geodesics);
  }

  for (auto &thread : threads) {
    if (thread.joinable())
      thread.join();
  }
}

template <typename _FullGeodesicData, typename _Spacetime, typename DiskTex>
RGBd get_single_geodesic_color(
    const _Spacetime &spacetime,
#if SKY_ENABLED
    const Image &sky,
#endif
    const DiskTex &disk_tex,
    const _FullGeodesicData &full) {
  (void)disk_tex;
  (void)spacetime;

  const int dead_reason = full.extra.dead_reason;

  if (dead_reason == DEAD_BLACK_HOLE) {
    return BLACK_HOLE_COLOR;
  } else if (dead_reason == DEAD_OBJECT) {
    return OBJECT_COLOR;
#if RENDER_DISK
  } else if (dead_reason == DEAD_DISK) {

    return disk_tex.get_color(disk_tex.get_tex_coord(spacetime, full));
#endif
#if FAKE_SKY
  } else if (dead_reason >= DEAD_SKY_TEX_OFFSET) {
    typedef typename _FullGeodesicData::basic_type basic_type;
    typedef typename basic_type::coord_type coord_type;
    auto spherical = full.basic.position.spherical_part(
          spacetime.coord_system_parameters(coord_type()));
    spherical.phi += M_PI;
    return RGBd{
      .5 + .5 * std::sin(spherical.theta) * std::cos(spherical.phi),
      .5 + .5 * std::sin(spherical.theta) * std::sin(spherical.phi),
      .5 + .5 * std::cos(spherical.theta)
    };
#endif
  } else {
#if SKY_ENABLED
    std::tie(tx, ty) = basic.position.spherical_part(
          spacetime.coord_system_parameters(_Coord())).to_txty();
    return sky.get_pixel_rel(tx, ty).to_RGBd();
#else
    return SKY_COLOR;
#endif
  }
}


template <typename _Spacetime, typename _FullGeodesicData, typename DiskTex>
void colorize_from_matrix_snapshot(
    const SnapshotMatrix<_FullGeodesicData> &snapshot,
    const _Spacetime &spacetime,
#if SKY_ENABLED
    const Image &sky,
#endif
    const DiskTex &disk_tex,
    RGBd *output) {
  auto geodesics = snapshot.geodesics;

  for (int k = 0; k < snapshot.width * snapshot.height; ++k) {
    const auto &geo = geodesics[k];
    output[k] = get_single_geodesic_color(
        spacetime,
#if SKY_ENABLED
        sky,
#endif
        disk_tex,
        geo);
  }

#if RENDER_DISK == DISK_SHAKURA
  {
    double total = 0;
    double total_dnu = 0;
    for (int i = 0; i < snapshot.height; ++i) {
      for (int j = 0; j < snapshot.width; ++j) {
        const auto &geo = geodesics[i * snapshot.width + j];
        if (geo.extra.dead_reason == DEAD_DISK) {
          auto tex_coord = disk_tex.get_tex_coord(spacetime, geo);
          total += tex_coord.intensity;
          total_dnu += tex_coord.intensity_dnu;
        }
      }
    }
    fprintf(stderr, "    Total power: %lg W/m^2 * px_area\n", total / UNIT_W);
    fprintf(stderr, "Total power dnu: %lg dnu W/m^2 * px_area\n", total_dnu / UNIT_W);
    fprintf(stderr, "Luminosity: %lg W\n",
        disk_tex.calc_total_luminosity() / UNIT_W);
  }
#endif

#if DISK_POLARIZATION
  int total_steps_disk = 0;
  int total_steps_other = 0;
  int width = snapshot.width;
  int height = snapshot.height;
  // int step = (int)(0.31 * std::sqrt(width + height));
  // double arrow_length = 1.5 * step;
  int step = (int)(.7 * std::sqrt(width + height));
  double arrow_length = .2 * step;
  double total_Q = 0;
  double total_U = 0;
  double total_absQ = 0;
  double total_absU = 0;
  double total_I = 0;
  double min_g = 1e10, max_g = -1e10;
  int disk_area = 0;
  for (int i = 0; i < snapshot.height; ++i) {
    for (int j = 0; j < snapshot.width; ++j) {
      const auto &geo = geodesics[i * snapshot.width + j];
      if (geo.extra.dead_reason == DEAD_DISK) {
        ++disk_area;
        auto tex_coord = disk_tex.get_tex_coord(spacetime, geo);
        double dQ = tex_coord.delta * tex_coord.intensity * tex_coord.chi2_cos;
        double dU = tex_coord.delta * tex_coord.intensity * tex_coord.chi2_sin;
        total_I += tex_coord.intensity;

        total_Q += dQ;
        total_U += dU;
        total_absQ += std::abs(dQ);
        total_absU += std::abs(dU);
        min_g = std::min(min_g, tex_coord.redshift_inf);
        max_g = std::max(max_g, tex_coord.redshift_inf);

        if (i % step == 0 && j % step == 0) {
          auto arrow = disk_tex.get_arrow_vector(tex_coord);
          // render_arrow((double)j, (double)i,
          //              // 300 * dQ,
          //              // 300 * dU,
          //              arrow_length * arrow.first,
          //              arrow_length * arrow.second,
          //              arrow_length, 0.3,
          //              ARROW_COLOR,
          //              width, height, output);
          double dx = arrow_length * arrow.first;
          double dy = arrow_length * arrow.second;
          render_line(j - dx, i - dy, j + dx, i + dy,
                      ARROW_COLOR, width, height, output);
        }
        total_steps_disk += geo.extra.steps;
      } else {
        total_steps_other += geo.extra.steps;
      }
    }
  }

  double total_chi = .5 * std::atan2(total_U, total_Q);
  double total_delta = std::abs(total_U / (total_I * std::sin(2 * total_chi)));
  fprintf(stderr, "Disk results:\n");
  fprintf(stderr, "  Debug min max g: %lg %lg\n", min_g, max_g);
  fprintf(stderr, "  Total absQ: %lg\n", total_absQ);
  fprintf(stderr, "  Total absU: %lg\n", total_absU);
  fprintf(stderr, "  Total Q: %lg\n", total_Q);
  fprintf(stderr, "  Total U: %lg\n", total_U);
  fprintf(stderr, "  Total I: %lg\n", total_I);
  fprintf(stderr, "  I/px: %lg\n", total_I / disk_area);
  fprintf(stderr, "  Total chi: %lg = %lg deg\n",
      total_chi, total_chi * 180 / M_PI);
  fprintf(stderr, "  Total delta: %lg\n", total_delta);

  double _scale = disk_tex._get_arrow_scale_01();
  double dx = arrow_length * _scale;
  double dy = 0;
  double x = .9 * width;
  double y = .9 * height;
  render_line(x, y, x + 2 * dx, y + 2 * dy, ARROW_COLOR, width, height, output);

  // render_arrow(.2 * width, .2 * height,
  //              total_Q,
  //              total_U,
  //              arrow_length, 0.3,
  //              ARROW_COLOR,
  //              snapshot.width, snapshot.height, output);


  fprintf(stderr, "Steps disk: %d  other: %d --> useful=%.2lf%%\n",
      total_steps_disk, total_steps_other,
      100. * total_steps_disk / (total_steps_disk + total_steps_other));
#endif

}


#endif
