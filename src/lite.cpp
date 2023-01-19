// Lite version of the code.
//
// Handles only the geodesic integration, with no polarization nor QED effects.

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <tuple>

enum Target {
  TARGET_BLACK_HOLE,
  TARGET_DISK,
  TARGET_INFINITY,
  TARGET_TOO_MANY_STEPS
};

/* For all physical constants, we assume G = c = 1! */
// ORIGINAL: include/physical_constants.h, include/parameters.h
constexpr double M = 0.05;
constexpr double a = 0.998 * M;
constexpr double aa = a * a;
constexpr double rs = 2 * M;  // rs = 2 G M / c^2

constexpr double disk_inner_r = 3 * rs;
constexpr double disk_outer_r = 10 * rs;

constexpr int MAX_N = 100000;
constexpr double epsilon = 1e-6;
constexpr double dlambda_initial = rs;
constexpr double dlambda_min = rs / 100000;
constexpr double dlambda_max = rs * 100;
constexpr double max_r = 1100 * rs;

/* Helper functions. */
template <typename T> inline T sqr(const T &x) { return x * x; }


/******************************** STRUCTS ***********************************/

struct BoyerLindquist {  // ORIGINAL: include/coordinate.h
  double t;
  double r;
  double theta;
  double phi;

  inline double &operator[](int k) { return ((double *)this)[k]; }
  inline double operator[](int k) const { return ((const double *)this)[k]; }

  /* Arithmetic operations on vectors, required for integration. */
  friend inline BoyerLindquist operator+(const BoyerLindquist &A,
                                         const BoyerLindquist &B) {
    return {A.t + B.t, A.r + B.r, A.theta + B.theta, A.phi + B.phi};
  }
  friend inline BoyerLindquist operator*(double c, const BoyerLindquist &B) {
    return {c * B.t, c * B.r, c * B.theta, c * B.phi};
  }
  inline BoyerLindquist operator-() const {
    return {-t, -r, -theta, -phi};
  }

  /* Estimate of the vector difference (error). */
  friend inline double numerical_distance(const BoyerLindquist &A,
                                          const BoyerLindquist &B) {
    return std::sqrt(sqr(A.t - B.t)
                   + sqr(A.r - B.r)
                   + sqr(A.theta - B.theta)
                   + sqr(A.phi - B.phi));
  }

  inline bool isfinite() const {
    using std::isfinite;
    return isfinite(t) && isfinite(r) && isfinite(theta) && isfinite(phi);
  }
};

struct Vector3 {  // ORIGINAL: include/vector.h
  double x, y, z;

  friend inline Vector3 operator+(const Vector3 &A, const Vector3 &B) {
    return {A.x + B.x, A.y + B.y, A.z + B.z};
  }
  friend inline Vector3 operator-(const Vector3 &A, const Vector3 &B) {
    return {A.x - B.x, A.y - B.y, A.z - B.z};
  }
  friend inline Vector3 operator*(double c, const Vector3 &B) {
    return {c * B.x, c * B.y, c * B.z};
  }
  inline Vector3 operator-() const {
    return {-x, -y, -z};
  }
  inline Vector3 normalized() const {
    return (1 / std::sqrt(x * x + y * y + z * z)) * (*this);
  }

  friend inline Vector3 cross(const Vector3 &A, const Vector3 &B) {
    return {
        A.y * B.z - A.z * B.y,
        A.z * B.x - A.x * B.z,
        A.x * B.y - A.y * B.x
    };
  }
};

struct Christoffel {
  double chr[4][4][4];
};

struct BGRA {
  unsigned char b, g, r, a;
};


/****************************** KERR METRIC **********************************/

Christoffel get_christoffel_ull(const BoyerLindquist &position) {
  /* Thomas Müller, Catalogue of Spacetimes, http://arxiv.org/abs/0904.4184v3 */
  /* ORIGINAL: include/kerr.h */
  double r = position.r;
  double rr = sqr(r);
  double cos_theta = std::cos(position.theta);
  double sin_theta = std::sin(position.theta);
  double cot_theta = cos_theta / sin_theta;
  double sin_cos_theta = sin_theta * cos_theta;
  double cos2_theta = sqr(cos_theta);
  double sin2_theta = sqr(sin_theta);
  double sigma = rr - aa * cos2_theta;
  double rr_aa = rr + aa;
  double one_over_Sigma = 1 / (rr + aa * cos2_theta);
  double one_over_Sigma2 = sqr(one_over_Sigma);
  double one_over_Sigma3 = one_over_Sigma * one_over_Sigma2;
  double Delta = rr - rs * r + aa;
  double one_over_Delta = 1 / Delta;

  /* Notation: t = t, r = r, h = theta, p = phi. */
  double rtt = .5 * rs * Delta * sigma * one_over_Sigma3;
  double ttr = .5 * rs * rr_aa * sigma * one_over_Sigma2 * one_over_Delta;
  double tth = -rs * aa * r * sin_cos_theta * one_over_Sigma2;
  double rtp = -.5 * rs * Delta * a * sin2_theta * sigma * one_over_Sigma3;
  double rrr = (r * aa * sin2_theta - .5 * rs * sigma) * one_over_Sigma
      * one_over_Delta;
  double rrh = - aa * sin_cos_theta * one_over_Sigma;
  double rhh = - r * Delta * one_over_Sigma;
  double php = cot_theta * (1 + rs * aa * r * sin2_theta * one_over_Sigma2);
  double htt = -rs * aa * r * sin_cos_theta * one_over_Sigma3;
  double ptr = .5 * rs * a * sigma * one_over_Sigma2 * one_over_Delta;
  double pth = -rs * a * r * cot_theta * one_over_Sigma2;
  double htp = rs * a * r * rr_aa * sin_cos_theta * one_over_Sigma3;
  double hrr = aa * sin_cos_theta * one_over_Sigma * one_over_Delta;
  double hrh = r * one_over_Sigma;
  double hhh = -aa * sin_cos_theta * one_over_Sigma;
  double thp = rs * a * aa * r * sin2_theta * sin_cos_theta * one_over_Sigma2;
  double trp = .5 * rs * a * sin2_theta * (aa * cos2_theta * (aa - rr)
      - rr * (aa + 3 * rr)) * one_over_Sigma2 * one_over_Delta;
  double prp = (r + .5 * rs * ((sqr(aa * sin_cos_theta) - rr * rr_aa)
      * one_over_Sigma2 - rr * one_over_Sigma)) * one_over_Delta;
  double rpp = Delta * sin2_theta * (-r
      + .5 * rs * aa * sin2_theta * sigma * one_over_Sigma2) * one_over_Sigma;
  double hpp = -sin_cos_theta * (sqr(rr_aa) - aa * Delta * sin2_theta
      + rr_aa * rs * aa * r * sin2_theta * one_over_Sigma) * one_over_Sigma2;

  return Christoffel{{
    {
      {0, ttr, tth, 0},
      {ttr, 0, 0, trp},
      {tth, 0, 0, 0},
      {0, trp, thp, 0},
    }, {
      {rtt, 0, 0, rtp},
      {0, rrr, rrh, 0},
      {0, rrh, rhh, 0},
      {rtp, 0, 0, rpp},
    }, {
      {htt, 0, 0, htp},
      {0, hrr, hrh, 0},
      {0, hrh, hhh, 0},
      {htp, 0, 0, hpp},
    }, {
      {0, ptr, pth, 0},
      {ptr, 0, 0, prp},
      {pth, 0, 0, php},
      {0, prp, php, 0},
    }
  }};
}


/**************************** GEODESIC EQUATION ******************************/

struct ODEState {  /* The 8D vector entering the diff. eq. */
  BoyerLindquist position;
  BoyerLindquist direction;

  friend inline ODEState operator+(const ODEState &A, const ODEState &B) {
    return {A.position + B.position, A.direction + B.direction};
  }
  friend inline ODEState operator*(double c, const ODEState &B) {
    return {c * B.position, c * B.direction};
  }
};

/* Compute y_{n+1} given y_n, fixed h and the RHS f(y). */
template <typename Vector, typename RHSFunc>
std::pair<Vector, Vector> integration_step_RKF45(
    const RHSFunc &RHS, const double h, const Vector &u) {
  // ORIGINAL: include/integration.h
  Vector k1, k2, k3, k4, k5, k6, out4, out5;

  k1 = h * RHS(u);
  k2 = h * RHS(u + (1. / 4) * k1);
  k3 = h * RHS(u + (3. / 32) * k1 + (9. / 32) * k2);
  k4 = h * RHS(u + (1932. / 2197) * k1
                 + (-7200. / 2197) * k2
                 + (7296. / 2197) * k3);
  k5 = h * RHS(u + (439. / 216) * k1
                 + (-8) * k2
                 + (3680. / 513) * k3
                 + (-845. / 4104) * k4);
  k6 = h * RHS(u + (-8. / 27) * k1
                 + 2 * k2
                 + (-3544. / 2565) * k3
                 + (1859. / 4104) * k4
                 + (-11. / 40) * k5);

  out4 = u + (25. / 216) * k1
           + (1408. / 2565) * k3
           + (2197. / 4104) * k4
           + (-1. / 5) * k5;
  out5 = u + (16. / 135) * k1
           + (6656. / 12825) * k3
           + (28561. / 56430) * k4
           + (-9. / 50) * k5
           + (2. / 55) * k6;

  return std::make_pair(out4, out5);  // Result of the 4th and 5th order.
}

/* The RHS of the diff. eq. dy/dlambda = f(y). Note: for simplicity, in the
 * code we use dlambda > 0 so we flipped here the signs for f(y). */
ODEState geodesic_RHS(const ODEState &state) {  // ORIGINAL: include/geodesic.h
  Christoffel christoffel_ull = get_christoffel_ull(state.position);

  ODEState result;
  for (int k = 0; k < 4; ++k) {
    double outer = 0;
    for (int i = 0; i < 4; ++i) {
      double inner = 0;
      for (int j = 0; j < 4; ++j)
        inner += christoffel_ull.chr[k][i][j] * state.direction[j];
      outer += state.direction[i] * inner;
    }
    result.direction[k] = outer;
  }
  result.position = -state.direction;
  return result;
}

/* One step of RKF45 integration, including h adaptivity. */
std::pair<ODEState, double> advance_geodesic_RKF45(const double min_h,
                                                   double h,
                                                   const double max_h,
                                                   const double epsilon,
                                                   const ODEState &state) {
  // ORIGINAL: include/raytracer.h
  // Here the h adaptivity is somewhat modified compared to RKF45.
  constexpr double SAFETY = 0.84;
  constexpr double SAFETY_INC = 1.2;
  constexpr double SAFETY_INC_MAX = 3.0;

  int limit = 0;
  for (;;) {
    ODEState final4, final5;
    std::tie(final4, final5) = integration_step_RKF45(geodesic_RHS, h, state);
    double R = (numerical_distance(final4.position, final5.position)
              + numerical_distance(final4.direction, final5.direction)) / h;
    double delta = SAFETY * std::pow(epsilon / R, 0.25);
    if (R <= epsilon) {
      double new_h = delta > SAFETY_INC_MAX
          ? std::min(h * SAFETY_INC_MAX, max_h)
          : (delta > SAFETY_INC ? std::min(h * delta / SAFETY_INC, max_h) : h);
      return std::make_pair(final4, new_h);
    } else if (h == min_h) {
      return std::make_pair(final4, h);
    }

    if (++limit == 20) {
      fprintf(stderr, "Too large step count in advance_geodesic_RKF45!\n");
      exit(1);
    }

    h *= delta;
    if (h < min_h) h = min_h;
    if (h > max_h) h = max_h;
  }
}

/* Integrates dy/dlambda = f(y), with all termination criteria. */
std::pair<ODEState, int> generate_geodesic(const BoyerLindquist &position,
                                           const BoyerLindquist &direction) {
  // ORIGINAL: include/raytracer.h
  ODEState state0{position, direction};  // Current y(lambda).
  double dlambda0 = dlambda_initial;     // The step dlambda.

  for (int n = 0; n < MAX_N; ++n) {
    ODEState state1;  // The state and the new stap dlambda.
    double dlambda1;
    std::tie(state1, dlambda1) = advance_geodesic_RKF45(
        dlambda_min, dlambda0, dlambda_max, epsilon, state0);

    if (state1.position.r > max_r)
      return std::make_pair(state1, TARGET_INFINITY);
    if (state1.position.r < rs * (1 + 1e-5)
        || !state1.position.isfinite()
        || !state1.direction.isfinite()) {
      return std::make_pair(state1, TARGET_BLACK_HOLE);
    }

    // Slowdown near the disk if we are crossing the z == 0 plane.
    const auto z0 = std::cos(state0.position.theta) * state0.position.r;
    const auto z1 = std::cos(state1.position.theta) * state1.position.r;
    const auto r = state0.position.r;
    if (((z0 >= 0 && z1 < 0) || (z0 <= 0 && z1 > 0))
        && 0.5 * disk_inner_r < r && r < 1.5 * disk_outer_r) {
      if (dlambda0 == dlambda_min) {
        if (disk_inner_r < r && r < disk_outer_r)
          return std::make_pair(state1, TARGET_DISK);
      } else if (0.7 * disk_inner_r < r && r < 1.3 * disk_outer_r) {
        dlambda0 = std::max(dlambda_min, dlambda0 / 2);
        continue;  // Decrease dlambda and repeat the integration step.
      }
    }

    state0 = state1;
    dlambda0 = dlambda1;
  }

  return std::make_pair(state0, TARGET_TOO_MANY_STEPS);
}


/********************************* CAMERA ************************************/

/* Computes the ray directions given rel. positions x, y (-1 <= x, y <= 1). */
class ProjectionCamera {  // ORIGINAL: include/raytracer.h
 private:
  Vector3 view, right;
  double fovy_tan;
 public:
  Vector3 eye, center, up;
  double aspect;

  ProjectionCamera(const Vector3 &_eye, const Vector3 &_center,
                   const Vector3 &_up, double _aspect, double _fovy)
        : eye(_eye), center(_center), up(_up), aspect(_aspect) {
    view = (center - eye).normalized();
    right = cross(view, up).normalized();
    up = cross(right, view);
    fovy_tan = std::tan(_fovy / 180 * M_PI / 2.0);
  }

  std::pair<Vector3, Vector3> ray(double x, double y) const {
    Vector3 direction = view
                      + (fovy_tan * aspect * x) * right
                      - (fovy_tan * y) * up;
    return std::make_pair(eye, direction);
  }
};

void matrix3_inverse(const double m[3][3], double out[3][3]) {
  // ORIGINAL: include/matrix.h
  double det = m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2])
             - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
             + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
  double invdet = 1 / det;
  out[0][0] = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * invdet;
  out[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invdet;
  out[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invdet;
  out[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invdet;
  out[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invdet;
  out[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * invdet;
  out[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invdet;
  out[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * invdet;
  out[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invdet;
}

/* Convert position and direction from Certesian to Boyer-Lindquist coord. */
std::pair<BoyerLindquist, BoyerLindquist> cartesian_to_boyer_lindquist(
    const Vector3 &pos, const Vector3 &dir) {  // ORIGINAL: include/coordinate.h
  double ss = sqr(pos.x) + sqr(pos.y);
  double s = std::sqrt(ss);
  double zz = sqr(pos.z);
  double RR_aa = ss + zz - aa;
  double rr = (RR_aa + std::sqrt(sqr(RR_aa) + 4 * aa * zz)) / 2;
  double r = std::sqrt(rr);
  double ra = std::sqrt(rr + aa);
  double cost = pos.z / r;
  double sint = s / ra;
  double cosp = pos.x / s;
  double sinp = pos.y / s;
  double inv[3][3], mat[3][3] = {
    {r / ra * sint * cosp, r * cost * cosp, -r * sint * sinp},
    {r / ra * sint * sinp, r * cost * sinp, r * sint * cosp},
    {cost, -r * sint, 0}
  };
  matrix3_inverse(mat, inv);

  BoyerLindquist position, direction;
  position.t = 0.0;  // Ignored.
  position.r = r;
  position.theta = std::atan2(s * r, ra * pos.z);
  position.phi = std::atan2(pos.y, pos.x);

  direction.t = 1.0;  // Approx., ok if the camera is far from the black hole.
  direction.r     = inv[0][0] * dir.x + inv[0][1] * dir.y + inv[0][2] * dir.z;
  direction.theta = inv[1][0] * dir.x + inv[1][1] * dir.y + inv[1][2] * dir.z;
  direction.phi   = inv[2][0] * dir.x + inv[2][1] * dir.y + inv[2][2] * dir.z;

  return std::make_pair(position, direction);
}


/********************* GENERATING AND SAVING THE IMAGE ***********************/

/* TGA (in the basic settings) is quite simple to handle.
 * Note: colors are stored in the BGRA order, not RGBA! */
void save_TGA(const BGRA *image, int width, int height, const char *filename) {
  // ORIGINAL: src/tga.cpp
  unsigned char header[18] = {0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  header[12] = width         & 0xFF;
  header[13] = (width >> 8)  & 0xFF;
  header[14] = height        & 0xFF;
  header[15] = (height >> 8) & 0xFF;
  header[16] = 32;  // Bits per pixel.
  header[17] = 32;  // (Disable) vertical mirroring.

  FILE *f = fopen(filename, "wb");
  if (f == nullptr) exit(10);
  if (fwrite(header, 18, 1, f) != 1) exit(11);
  if (fwrite(image, width * height * 4, 1, f) != 1) exit(12);
  fclose(f);
}

/* Pixel color, depending on the target hit and the final geodesic state. */
BGRA get_color(const ODEState &state, int target) {
  // ORIGINAL: include/render.h
  switch(target) {
    case TARGET_BLACK_HOLE:
      return BGRA{0, 0, 255, 255};
    case TARGET_DISK: {
      double phi = 24 * state.position.phi / (2 * M_PI);
      phi = 2 * M_PI * std::floor(phi) / 24;
      double radial = (disk_outer_r - state.position.r)
                    / (disk_outer_r - disk_inner_r);
      unsigned char r = 255 * (.5 - .5 * std::cos(phi));
      unsigned char g = 255 * (.25 + (int)(radial * 10) / 10. * .75);
      return BGRA{0, g, r, 255};
    }
    case TARGET_INFINITY:
      return BGRA{0, 0, 0, 255};
    // case TARGET_TOO_MANY_STEPS:
    default:
      return BGRA{255, 0, 255, 255};
  };
}

int main() {
  const int width = 640;
  const int height = 3 * width / 4;
  auto image = std::make_unique<BGRA[]>(width * height);
  if (image == nullptr) return 1;

  double fovy = 1.0;                       // Field of view.
  double camera_dist = 1000 * rs;          // Distance from the origin.
  double camera_theta = 75 * M_PI / 180.;  // Angle of the position wrt z-axis.
  double cos_theta = std::cos(camera_theta);
  double sin_theta = std::sin(camera_theta);
  assert(camera_dist < max_r);
  // Camera's position (eye), its orientation (up) and the target (center).
  Vector3 eye = camera_dist * Vector3{sin_theta, 0, cos_theta};
  Vector3 up = camera_dist * Vector3{-cos_theta, 0, sin_theta};
  Vector3 center{0, 0, 0};
  ProjectionCamera camera(eye, center, up, (double)width / height, fovy);

  // ORIGINAL: include/render.h (generate_image)
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      double x = double(1 + 2 * j - width) / width;    // -1 <= x <= 1
      double y = double(1 + 2 * i - height) / height;  // -1 <= y <= 1
      std::pair<Vector3, Vector3> ray = camera.ray(x, y);

      BoyerLindquist position, direction;
      std::tie(position, direction) = cartesian_to_boyer_lindquist(
          ray.first, -ray.second);  // Minus because we send the ray backwards.

      std::pair<ODEState, int> result = generate_geodesic(position, direction);
      image[i * width + j] = get_color(result.first, result.second);
    }
    if (i % 16 == 0) fprintf(stderr, "%.1lf%% ", 100. * (i + 1) / height);
  }
  fprintf(stderr, "\n");

  save_TGA(image.get(), width, height, "output_lite.tga");
  return 0;
}
