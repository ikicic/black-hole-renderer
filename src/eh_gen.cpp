#include "base.h"
#include "euler_heisenberg.h"
#include "qed_lagrangian.h"
// #include "flat.h"  // Disabled.
#include "raytracer.h"
#include "autodiff.h"

#include "3rd/sha1.h"

#include <chrono>
#include <thread>
#include <sstream>

namespace QEDCache {

/* Everything is in the units of B_c. */
constexpr int CNT_F = 1000000;
constexpr int CNT_G = 1;
constexpr double MIN_F = sqr(1e-11);
constexpr double MIN_G = sqr(1e-11);
constexpr double MAX_F = sqr(300);
constexpr double MAX_G = sqr(1e-11);
constexpr int SIZE = CNT_F * CNT_G;

LambdaCacheType *lambda_cache;

#if 0
std::pair<double, double> coord_to_FG(int j, int i) {
  return std::make_pair(
      MIN_F * std::exp(std::log(MAX_F / MIN_F) * (double)j / (CNT_F - 1)),
      MIN_G * std::exp(std::log(MAX_G / MIN_G) * (double)i / (CNT_G - 1))
  );
}
std::pair<double, double> FG_to_coord(double F, double G) {
  return std::make_pair(
      std::log(F / MIN_F) / std::log(MAX_F / MIN_F) * (CNT_F - 1),
      std::log(G / MIN_G) / std::log(MAX_G / MIN_G) * (CNT_G - 1)
  );
}

LambdaCacheType get_cached_lambdas__dimless(double F, double G) {
  double x, y;
  std::tie(x, y) = FG_to_coord(F + 1e-20, G + 1e-20);
  int i = (int)y;
  int j = (int)x;
  //if (i < 0 || j < 0 || i >= CNT_G - 1 || j >= CNT_F - 1) {

  // Anyway inside the neutron star...
  if (i >= CNT_G - 1)
    i = CNT_G - 1;
  if (j >= CNT_F - 1)
    j = CNT_F - 1;
  if (i < 0)
    i = 0;

  if (i < 0 || j < 0) {
    fprintf(stderr, "Cache boundaries exceeded.\n");
    fprintf(stderr, "i=%d j=%d\n", i, j);
    fprintf(stderr, "MAX_F=%lg  F=%lg  MIN_F=%lg\n", MAX_F, F, MIN_F);
    fprintf(stderr, "MAX_G=%lg  G=%lg  MIN_G=%lg\n", MAX_G, G, MIN_G);
    fprintf(stderr, "B=%lg T\n", std::sqrt(F) * PHY_Bc / UNIT_T);
    exit(1);
  }


  // std::cerr << "\n";
  // std::cerr << "i=" << i << "  j=" << j << '\n';
  // std::cerr << "x=" << x << "  y=" << y << '\n';
  // std::cerr << lambda_cache[i * CNT_F + j].lambda[1] << '\n';
  // std::cerr << lambda_cache[i * CNT_F + j + 1].lambda[1] << '\n';
  // std::cerr << lambda_cache[(i + 1) * CNT_F + j + 1].lambda[1] << '\n';
  // std::cerr << lambda_cache[(i + 1) * CNT_F + j].lambda[1] << '\n';

  x -= j;
  y -= i;

  return (1 - x) * (1 - y) * lambda_cache[i * CNT_F + j]
       + x * (1 - y) * lambda_cache[i * CNT_F + j + 1]
       + x * y * lambda_cache[(i + 1) * CNT_F + j + 1]
       + (1 - x) * y * lambda_cache[(i + 1) * CNT_F + j];
}

#endif

inline double coord_to_F(int j) {
  return MIN_F * std::exp(std::log(MAX_F / MIN_F) * (double)j / (CNT_F - 1));
}
inline double F_to_coord(double F) {
  return std::log(F / MIN_F) / std::log(MAX_F / MIN_F) * (CNT_F - 1);
}

LambdaCacheType get_cached_lambdas__dimless(double F, double /* G */) {
  double x = F_to_coord(F + 1e-20);
  int j = (int)x;

  // Anyway inside the neutron star...
  if (j >= CNT_F - 1)
    j = CNT_F - 1;

  if (j < 0) {
    fprintf(stderr, "Cache boundaries exceeded.\n");
    fprintf(stderr, "j=%d\n", j);
    fprintf(stderr, "MAX_F=%lg  F=%lg  MIN_F=%lg\n", MAX_F, F, MIN_F);
    fprintf(stderr, "B=%lg T\n", std::sqrt(F) * PHY_Bc / UNIT_T);
    exit(1);
  }

  x -= j;
  return (1 - x) * lambda_cache[j] + x * lambda_cache[j + 1];
}

LambdaCacheType get_cached_lambdas(double F, double G) {
  LambdaCacheType result = get_cached_lambdas__dimless(
      F * PHY_inv_sqr_Bc,
      G * PHY_inv_sqr_Bc);
  for (int i = 0; i < 2; ++i) {
    result.lambda[i].value() *= PHY_inv_sqr_Bc;
    result.lambda[i].first(0) *= sqr(PHY_inv_sqr_Bc);
    result.lambda[i].first(1) *= sqr(PHY_inv_sqr_Bc);
  }
  return result;
}

void _lambda_preprocess_status_thread(std::vector<int> &status, int total) {
  auto start_time = std::chrono::system_clock::now();
  int done;
  do {
    done = 0;
    for (int x : status)
      done += x;
    std::chrono::duration<double> time_delta =
        std::chrono::system_clock::now() - start_time;
    int ETA = (int)((total - done) * time_delta.count() / done);
    fprintf(stderr, "%.2lf%% (%.0lfOPS, ETA: %ds)      \r",
        100. * done / total,
        done / time_delta.count(),
        ETA
      );
    std::this_thread::sleep_for(std::chrono::milliseconds(
          ETA < 2 ? 100 : 500));
  } while (done < total);
  std::chrono::duration<double> time_delta =
      std::chrono::system_clock::now() - start_time;
  fprintf(stderr, "\nTotal time: %.1lfs\n", time_delta.count());
}

LambdaCacheType _calc_lambdas__dimless(double F) {
  lambda_t Fad(F, 1, 0);
  lambda_t Gad(MIN_G, 0, 1);

  std::pair<lambda_t, lambda_t> result = qed_metric_correction_lambda__dimless(
      [](auto F, auto G) { return EH::lagrangian_real__dimless(F, G); },
      Fad, Gad);
  // std::cerr << result.first << '\n';
  // std::cerr << result.second << '\n';
  return {result.first, result.second};
}

void _lambda_preprocess_worker_thread(int thread_cnt, int index, int *status) {
  *status = 0;
  for (int k = index; k < SIZE; k += thread_cnt) {
#if 0
    int i = k / CNT_F;
    int j = k % CNT_F;
    double F, G;
    std::tie(F, G) = coord_to_FG(j, i);
    lambda_t Fad(F, 1, 0);
    lambda_t Gad(G, 0, 1);
#endif
    int j = k;
    double F = coord_to_F(j);
    lambda_cache[k] = _calc_lambdas__dimless(F);
    ++*status;
  }
}

bool _lambda_preprocess(int thread_count) {
  if (lambda_cache != nullptr)
    return true;
  std::stringstream s;
  s << CNT_F << "|" << CNT_G << "|";
  s << MAX_F << "|" << MAX_G << "|";
  s << MIN_F << "|" << MIN_G << "|";
  std::string filename = "output/preprocess/" + sha1(s.str()) + ".dat";
  std::cerr << "CACHE: " << filename << '\n';

  lambda_cache = new LambdaCacheType[SIZE];
  FILE *f = fopen(filename.c_str(), "rb");
  if (f != nullptr) {
    int loaded = fread(lambda_cache, sizeof(LambdaCacheType), SIZE, f);
    fclose(f);
    if (loaded != SIZE) {
      fprintf(stderr, "Error loading lambda cache file! Aborting!\n");
      delete []lambda_cache;
      lambda_cache = nullptr;
      return false;
    }
    fprintf(stderr, "Lambda cache file loaded!\n");
    return true;
  }

  std::cerr << s.str() << '\n';
  fprintf(stderr, "Lambda cache file %s not found, generating.\n",
      filename.c_str());
  std::vector<std::thread> threads;
  std::vector<int> status(thread_count);
  for (int i = 0; i < thread_count; ++i) {
    threads.emplace_back(_lambda_preprocess_worker_thread,
                         thread_count,
                         i,
                         &status[i]);
  }
  threads.emplace_back(_lambda_preprocess_status_thread,
                       std::ref(status),
                       SIZE);

  for (auto &thread : threads)
    if (thread.joinable())
      thread.join();

  fprintf(stderr, "  Saving...");
  f = fopen(filename.c_str(), "wb");
  if (fwrite(lambda_cache, sizeof(LambdaCacheType), SIZE, f) != SIZE) {
    fprintf(stderr, "Error writing lambda cache! Aborting!\n");
    delete []lambda_cache;
    lambda_cache = nullptr;
    return false;
  }
  fclose(f);

  fprintf(stderr, "  Done!\n");
  return true;
}

bool lambda_preprocess(int thread_count) {
  bool result = _lambda_preprocess(thread_count);
  if (result && debug) {
    double worst = 0;
    for (int j = 1; j < CNT_F; ++j) {
      double ratio = lambda_cache[j].lambda[1].value()
                   / lambda_cache[j - 1].lambda[1].value();
      worst = std::max(std::abs(ratio - 1), worst);
    }
    std::cerr << "Worst rel difference: " << worst << '\n';

    for (int i = 0; i < std::min(3, CNT_G); ++i) {
      for (int j = 0; j < 6; ++j) {
        std::cerr << j << " "
                  << i << " "
#if 0
                  << coord_to_FG(j, i).first << " B_c^2  "
                  << coord_to_FG(j, i).second << " B_c^2  "
#endif
                  << coord_to_F(j) << " B_c^2  "
                  << lambda_cache[i * CNT_F + j].lambda[1] << " B_c^-2\n";
      }
      std::cerr << '\n';
    }

    // for (int i = 0; i < 100; ++i) {
    //   typedef first_partial_derivatives<double, 2> fpd;
    //   double B = i * 1e-9 * PHY_Bc;
    //   double F = sqr(B) / 2;
    //   double G = 0;
    //   LambdaCacheType cached = get_cached_lambdas(F, G);
    //   auto calced = qed_metric_correction_lambda(
    //       [](auto F, auto G) { return EH::lagrangian_real(F, G); },
    //       fpd(F, 1, 0), fpd(G, 0, 1));
    //   // auto approx = qed_metric_correction_lambda(
    //   //     [](auto F, auto G) { return QED::lagrangian_lowest_order(F, G); },
    //   //     fpd(F, 1, 0), fpd(G, 0, 1));
    //   cached.lambda[1].value() *= sqr(PHY_Bc);
    //   cached.lambda[1].first(0) *= sqr(sqr(PHY_Bc));
    //   cached.lambda[1].first(1) *= sqr(sqr(PHY_Bc));
    //   calced.second.value() *= sqr(PHY_Bc);
    //   calced.second.first(0) *= sqr(sqr(PHY_Bc));
    //   calced.second.first(1) *= sqr(sqr(PHY_Bc));
    //   std::cerr.precision(10);
    //   std::cerr << "B=" << B
    //             << "   B=" << B / UNIT_T << " T"
    //             << "\t" /* << cached.lambda[0] << " " */ << cached.lambda[1]
    //             << "\t" /* << calced.first     << " " */ << calced.second
    //             << '\n';
    //         //    << "\t" /* << approx.first     << " " */ << approx.second << '\n';
    // }
    // exit(1);
  }
  return result;
}

#if CHECK_LAMBDA_PRECISION
void check_lambda_interpolation_precision(void) {
  /* Check worst relative precision of interpolation on preprocessed data. */
  constexpr int N = 100;
  for (int k = 0; k < 100; ++k)
    std::cerr << coord_to_F(k) << "\t" << lambda_cache[k].lambda[1] << '\n';

  constexpr double min_f = 1e-6;
  LambdaCacheType worst{lambda_t(0), lambda_t(0)};
  for (int k = 0; k < N; ++k) {
    double F = min_f * std::pow(MAX_F / min_f, random_double(0, 1));
    // double F = min_f * std::pow(MAX_F / min_f, 0.000001 * k);
    double G = 0;
    LambdaCacheType val = get_cached_lambdas__dimless(F, G);
    LambdaCacheType exp = _calc_lambdas__dimless(F);

    auto update = [](auto &w, auto value, auto expected) {
      w = std::max(w, std::abs(value / expected - 1));
    };

    for (int i = 0; i < 2; ++i) {
      update(worst.lambda[i].value(), val.lambda[i].value(), exp.lambda[i].value());
      update(worst.lambda[i].first(0), val.lambda[i].first(0), exp.lambda[i].first(0));
      update(worst.lambda[i].first(1), val.lambda[i].first(1), exp.lambda[i].first(1));
    }
    // std::cerr << F << "\t" << val.lambda[1] << "\t" << exp.lambda[1] << "\t" << worst.lambda[1] << "\n";
    if (k % 10 == 9) fprintf(stderr, "%.2lf%% ", 100. * (k + 1) / N);
  }

  std::cerr << "\nWorst relative error: "
            << worst.lambda[0] << "  "
            << worst.lambda[1] << "\n";
}
#endif

};



namespace QED {


#if GENERATE_LAMBDAS
void generate_lambdas(void) {
#ifndef PREPROCESS_LAMBDAS
#error PREPROCESS_LAMBDAS must be enabled!
#endif
  // FILE *f = fopen("output/eh_data.csv", "w");
  FILE *f = fopen("output/eh_data_large.csv", "w");
  if (f == NULL) {
    fprintf(stderr, "Error opening gen output file.\n");
    exit(1);
  }
  constexpr int N = 100;
  for (int i = 0; i < N; ++i) {
    // double B = 70. * i / N;  // dimless.
    double B = 10000. * i / N;  // dimless.
    double F = sqr(B) / 2;
    // double G = 0;

    QEDCache::LambdaCacheType lambdas = QEDCache::_calc_lambdas__dimless(F);
    // QEDCache::LambdaCacheType lambdas = QEDCache::get_cached_lambdas__dimless(F, G);

    // fprintf(f, "%.10lg, %.10lg, %.10lg\n",
    //     B, lambdas.lambda[0].value(), lambdas.lambda[1].value());
    if (i % 10 == 9)
      fprintf(stderr, "%.2lf%% ", 100. * (i + 1) / N);
  }

  fclose(f);
  fprintf(stderr, "Gen done.\n");
}
#endif

};
