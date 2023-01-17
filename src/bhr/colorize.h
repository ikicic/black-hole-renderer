#ifndef COLORIZE_H
#define COLORIZE_H

#include <bhr/mod.h>

#define BLOOM  0

#include <bhr/render.h>
#include <bhr/recursive_render.h>


template <typename _FullGeodesicData, typename... Args>
void colorize(const Snapshot<_FullGeodesicData> *snapshot, Args... args) {
  // TODO: How to avoid this? Virtual templated functions are not supported.
  {
    auto p = dynamic_cast<const SnapshotMatrix<_FullGeodesicData> *>(snapshot);
    if (p) {
      colorize_from_matrix_snapshot(*p, args...);
      return;
    }
  }
  {
    auto p = dynamic_cast<const SnapshotRecursive<_FullGeodesicData> *>(snapshot);
    if (p) {
      colorize_from_recursive_snapshot(*p, args...);
      return;
    }
  }

  fprintf(stderr, "What is %p?\n", snapshot);
  exit(1);
}


/* TODO: Check these blurs, something seems not right. */

std::pair<double, double> gaussian_blur__get_mnu_and_factor(
    int radius, double intensity_at_edge) {
  double mnu = std::log(intensity_at_edge) / sqr(radius);  // Minus sigma.
  double total = 0;
  for (int x = -radius; x <= radius; ++x)
    total += std::exp(mnu * x * x);
  double global_factor = 1 / sqr(total);
  return std::make_pair(mnu, global_factor);
}

void gaussian_blur(RGBd *input, RGBd *output, RGBd *tmp, int width, int height,
    int radius, double intensity_at_edge) {
  double mnu, global_factor;
  std::tie(mnu, global_factor) = gaussian_blur__get_mnu_and_factor(
      radius, intensity_at_edge);

  for (int k = 0; k < width * height; ++k)
    tmp[k] = RGBd{0.0, 0.0, 0.0};
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      int k = i * width + j;
      for (int c = 0; c < 3; ++c) {
        double extra = input[k][c];
        if (!extra)
          continue;
        for (int y = std::max(-i, -radius); y <= radius; ++y) {
          if (i + y >= height)
            break;
          double factor = global_factor * std::exp(mnu * y * y);
          tmp[k + y * width][c] += factor * extra;
        }
      }
    }
  }

  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      int k = i * width + j;
      for (int c = 0; c < 3; ++c) {
        double extra = tmp[k][c];
        if (!extra)
          continue;
        for (int x = std::max(-j, -radius); x <= radius; ++x) {
          if (j + x >= width)
            break;
          double factor = std::exp(mnu * x * x);
          output[k + x][c] += factor * extra;
        }
      }
    }
  }
}

void selective_gaussian_blur(RGBd *input, double *sigma, RGBd *output,
    int width, int height) {
  /* Perform gaussian blur with the kernel equal to:
   * G(r) = (1 / sqrt(2 pi sigma^2)) exp(-x^2 / (2 sigma^2))
   */

  std::vector<double> erfs;

  for (int i = 0; i < height; ++i)
    for (int j = 0; j < width; ++j) {
      int k = i * width + j;
      if (sigma[k] == 0) {
        output[k] += input[k];
        continue;
      }

      double nu = std::sqrt(.5) / sigma[k];
      int radius = 1;
      while (nu * sqr(radius) < 300)
        ++radius;

      erfs.resize(2 * radius + 2);
      for (int l = -radius - 1; l <= radius; ++l)
        erfs[l + radius + 1] = .5 * std::erf(nu * (l + .5));

      for (int y = std::max(-i, -radius); y <= radius && i + y < height; ++y) {
        double factory = erfs[y + radius + 1] - erfs[y + radius];
        for (int x = std::max(-j, -radius); x <= radius && j + x < width; ++x) {
          double factorx = erfs[x + radius + 1] - erfs[x + radius];
          output[k + y * width + x] += (factorx * factory) * input[k];
        }
      }
    }
}

void _gaussian_blur2(RGBd *input, RGBd *output, RGBd *tmp,
    int width, int height, int from, int to,
    int radius, double mnu, double global_factor) {
  for (int i = from; i < to; ++i) {
    for (int j = 0; j < width; ++j) {
      int k = i * width + j;
      RGBd total{0, 0, 0};
      for (int y = std::max(-i, -radius); y <= radius && i + y < height; ++y)
        total += std::exp(mnu * y * y) * input[k + y * width];
      tmp[k] = global_factor * total;
    }
  }

  for (int i = from; i < to; ++i) {
    for (int j = 0; j < width; ++j) {
      int k = i * width + j;
      RGBd total{0, 0, 0};
      for (int x = std::max(-j, -radius); x <= radius && j + x < width; ++x)
        total += std::exp(mnu * x * x) * tmp[k + x];
      output[k] += total;
    }
  }

  // double ttmp = global_factor * std::exp(mnu);
  // for (int i = from; i < to; ++i)
  //   for (int j = 0; j < width; ++j) {
  //     int k = i * width + j;
  //     output[k] -= ttmp * input[k];
  //   }
}

inline void gaussian_blur2(RGBd *input, RGBd *output, RGBd *tmp,
    int width, int height, int radius, double intensity_at_edge) {
  double mnu, global_factor;
  std::tie(mnu, global_factor) = gaussian_blur__get_mnu_and_factor(
      radius, intensity_at_edge);

  _gaussian_blur2(
      input, output, tmp,
      width, height, 0, height, radius,
      mnu, global_factor);
}

void gaussian_blur__parallelized(RGBd *input, RGBd *output, RGBd *tmp,
    int width, int height, int radius, double intensity_at_edge,
    double additional_factor, int thread_count) {
  double mnu, global_factor;
  std::tie(mnu, global_factor) = gaussian_blur__get_mnu_and_factor(
      radius, intensity_at_edge);
  global_factor *= additional_factor;

  std::vector<std::thread> threads;
  for (int i = 0; i < thread_count; ++i) {
    int from = i * height / thread_count;
    int to = (i + 1) * height / thread_count;
    threads.emplace_back(
        _gaussian_blur2,
        input, output, tmp,
        width, height, from, to,
        radius, mnu, global_factor
    );
  }

  for (auto &thread : threads)
    if (thread.joinable())
      thread.join();
}

void apply_color_filters(
    RGBd *rgb, RGBA *rgba, int width, int height, int thread_count) {
  const int size = width * height;
  (void)size;
  (void)thread_count;

  // double max_I = 0;
  // for (int k = 0; k < size; ++k) {
  //   max_I = std::max(
  //       max_I,
  //       std::max(rgb[k][0], 0.) + std::max(rgb[k][1], 0.) + std::max(rgb[k][2], 0.)
  //   );
  // }

#if BLOOM
  RGBd *rgb_bright = new RGBd[size];
  RGBd *rgb_final = new RGBd[size];
  RGBd *rgb_tmp = new RGBd[size];

  double *sigma = new double[size];
  for (int k = 0; k < size; ++k) {
    rgb_bright[k] = rgb[k];
    sigma[k] = std::sqrt(std::log(1 + rgb[k].luminance()));
    rgb[k] = RGBd{0, 0, 0};
  }
  selective_gaussian_blur(rgb_bright, sigma, rgb, width, height);
  delete []sigma;

#define GAUSSIAN_BLUR(_radius, _intensity_at_edge, _factor) \
    gaussian_blur__parallelized( \
        rgb_bright, rgb, rgb_tmp, width, height, \
        (_radius), (_intensity_at_edge), _factor, thread_count);
  double a = .05;
  double b = .03;
  double c = .01;
  double factor = 1 / (1 + a + b + c);
  for (int k = 0; k < size; ++k)
    rgb[k] *= factor;
  GAUSSIAN_BLUR(width / 400, 1e-3, factor * a);
  GAUSSIAN_BLUR(width / 50, 1e-5, factor * b);
  GAUSSIAN_BLUR(width / 4, 1e-9, factor * c);
#undef GAUSSIAN_BLUR

  // ------- SELECTIVE GAUSSIAN BLUR TEST START ---------
  // double *sigma = new double[width * height];
  // for (int i = 0; i < height; ++i)
  //   for (int j = 0; j < width; ++j) {
  //     rgb[i * width + j] = RGBd{0, 0, 0};
  //     rgb_bright[i * width + j] = RGBd{(double)((i / 20 + (j + i) / 20) % 2), 0, 0};
  //     sigma[i * width + j] = 5 * (1 + std::sin(0.023 * (i + j * 0.1)));
  //   }
  // selective_gaussian_blur(rgb_bright, sigma, rgb, width, height);
  // -------- SELECTIVE GAUSSIAN BLUR TEST END ----------

  // for (int k = 0; k < size; ++k)
  //   rgb[k] += rgb_final[k];

  delete []rgb_tmp;
  delete []rgb_bright;
  delete []rgb_final;

  constexpr double GAMMA = .5;
  constexpr double A = 1;
#else
  constexpr double GAMMA = 1;
  constexpr double A = 1;
#endif

  constexpr int ZOOM = 1;
  for (int i = 0; i < height; ++i)
    for (int j = 0; j < width; ++j) {
      int k2 = (i / ZOOM) * width + j / ZOOM;
      // double factor = 255 * A * std::pow(rgb[k2].luminance(), GAMMA - 1);
      rgba[i * width + j] = RGBA(
          // clamp((int)(factor * rgb[k2][0]), 0, 255),
          // clamp((int)(factor * rgb[k2][1]), 0, 255),
          // clamp((int)(factor * rgb[k2][2]), 0, 255)
          clamp((int)(255 * A * std::pow(rgb[k2][0], GAMMA)), 0, 255),
          clamp((int)(255 * A * std::pow(rgb[k2][1], GAMMA)), 0, 255),
          clamp((int)(255 * A * std::pow(rgb[k2][2], GAMMA)), 0, 255)
      );
    }
}

#endif
