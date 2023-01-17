#ifndef BASE_H
#define BASE_H

#include <cmath>
#include <vector>

#ifndef __clang__
# define CMATH_CONSTEXPR  constexpr
#else
# define CMATH_CONSTEXPR
#endif

typedef double real_t;
typedef double colreal_t;  // color real

struct Null{ }; // empty class

extern int debug;

#define DISK_DUMMY      1
#define DISK_KERTAP     2
#define DISK_SHAKURA    3

#define CHECK_KERR          0
#define CHECK_LAMBDA_PRECISION  0
#define GENERATE_LAMBDAS    0

#define PREPROCESS_LAMBDAS  1
#define RENDER_DISK         DISK_KERTAP
#define DISK_POLARIZATION   (RENDER_DISK == DISK_KERTAP)
#define SKY_ENABLED         0
#define MAGNETIC_FIELD      0
#define MAGNETIC_FIELD_FULL 0

#if MAGNETIC_FIELD_FULL && !MAGNETIC_FIELD
#error MAGNETIC_FIELD_FULL requires MAGNETIC_FIELD
#endif

#define FAKE_SKY      0
#define DISK_RELIEF_TEXTURE 0

/* Params BEGIN. */
#define PREDEFINED_PARAMS 1

#if PREDEFINED_PARAMS
# define PARAMS_CONSTEXPR        constexpr
# define PARAMS_CMATH_CONSTEXPR  CMATH_CONSTEXPR
# define PARAMS_FUNC_CONST
# define PARAMS_FUNC_STATIC                  static
# define PARAMS_FUNC_STATIC_CONSTEXPR        static constexpr
# define PARAMS_FUNC_STATIC_CMATH_CONSTEXPR  static CMATH_CONSTEXPR
#else
/* blah... */
# define PARAMS_CONSTEXPR                    const
# define PARAMS_CMATH_CONSTEXPR              const
# define PARAMS_FUNC_CONST                   const
# define PARAMS_FUNC_STATIC                  inline
# define PARAMS_FUNC_STATIC_CONSTEXPR        inline
# define PARAMS_FUNC_STATIC_CMATH_CONSTEXPR  inline
#endif
/* Params END. */

#include <bhr/vector.h>
#include <bhr/coordinate.h>


#define DEAD_BLACK_HOLE     1
#define DEAD_FLAT           2
#define DEAD_DISK           3
#define DEAD_OBJECT         4
#define DEAD_UNUSED         30
#define DEAD_SKY_TEX_OFFSET 0x01000000



template <typename _T>
struct _RGB : Vector<_T, 3> {
  _RGB() {}
  _RGB(_T r, _T g, _T b) {
    (*this)[0] = r;
    (*this)[1] = g;
    (*this)[2] = b;
  }
  _RGB(const Vector<_T, 3> &parent) : Vector<_T, 3>(parent) {}

  template <typename _T2>
  explicit _RGB(const _RGB<_T2> &rgb) {
    (*this)[0] = rgb[0];
    (*this)[1] = rgb[1];
    (*this)[2] = rgb[2];
  }

  inline _T luminance(void) const {
    return _T(0.21) * (*this)[0]
         + _T(0.72) * (*this)[1]
         + _T(0.07) * (*this)[2];
  }
};

typedef Vector<colreal_t, 3> XYZd;
typedef _RGB<colreal_t> RGBd;
typedef _RGB<float> RGBf;

// utility.h
inline double numerical_sqr_distance(const RGBd &A, const RGBd &B) {
  return sqr(A[0] - B[0]) + sqr(A[1] - B[1]) + sqr(A[2] - B[2]);
}

struct RGBA {
  unsigned int rgba;

  RGBA() {}
  RGBA(unsigned int _rgba) : rgba(_rgba) {}
  RGBA(unsigned char r, unsigned char g, unsigned char b, unsigned char a = 255) {
    rgba = ((unsigned int)a << 24) | ((unsigned int)b << 16) | ((unsigned int)g << 8) | r;
  }

  inline int get_r(void) const { return rgba & 0xFF; }
  inline int get_g(void) const { return (rgba >> 8) & 0xFF; }
  inline int get_b(void) const { return (rgba >> 16) & 0xFF; }
  inline RGBA get_BGRA(void) const {
    return (rgba & 0xFF00FF00) | ((rgba & 0xFF) << 16) | ((rgba >> 16) & 0xFF);
  }
  inline RGBd to_RGBd(void) const {
    return RGBd(
      (1. / 255) * (rgba & 0xFF),
      (1. / 255) * ((rgba >> 8) & 0xFF),
      (1. / 255) * ((rgba >> 16) & 0xFF)
    );
  }
};


template <typename _T1, typename _T2, typename _T3>
inline bool is_between(const _T1 &x, const _T2 &low, const _T3 &high) {
  return low <= x && x <= high;
}
template <typename _T>
inline int int_sgn(const _T &x) {
  return x > 0 ? 1 : (x < 0 ? -1 : 0);
}
inline double random_double(double low, double high) {
  return low + (high - low) * rand() / RAND_MAX;
}

struct Image;

template <typename _FullGeodesicData>
class Snapshot {
 public:
  Snapshot(int _width, int _height) : width(_width), height(_height) {}
  virtual bool load(FILE *f) = 0;
  virtual bool save(FILE *f) const = 0;
  virtual void save_extra(void) const { };

  int width;
  int height;
};

#endif
