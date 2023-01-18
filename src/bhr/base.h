#ifndef BASE_H
#define BASE_H

#include <cmath>
#include <vector>

#include <bhr/config.h>
#include <bhr/vector.h>

namespace bhr {

extern int debug;

#define DEAD_BLACK_HOLE     1
#define DEAD_FLAT           2
#define DEAD_DISK           3
#define DEAD_OBJECT         4
#define DEAD_UNUSED         30
#define DEAD_SKY_TEX_OFFSET 0x01000000

template <typename T>
struct RGB : Vector<T, 3> {
  RGB() {}
  RGB(T r, T g, T b) {
    (*this)[0] = r;
    (*this)[1] = g;
    (*this)[2] = b;
  }
  RGB(const Vector<T, 3> &parent) : Vector<T, 3>(parent) {}

  template <typename T2>
  explicit RGB(const RGB<T2> &rgb) {
    (*this)[0] = rgb[0];
    (*this)[1] = rgb[1];
    (*this)[2] = rgb[2];
  }

  inline T luminance(void) const {
    return T(0.21) * (*this)[0]
         + T(0.72) * (*this)[1]
         + T(0.07) * (*this)[2];
  }
};

typedef Vector<colreal_t, 3> XYZd;
typedef RGB<colreal_t> RGBd;
typedef RGB<float> RGBf;

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


inline double random_double(double low, double high) {
  return low + (high - low) * rand() / RAND_MAX;
}

template <typename FullGeodesicData>
class Snapshot {
 public:
  Snapshot(int _width, int _height) : width(_width), height(_height) {}
  virtual ~Snapshot() = default;
  virtual bool load(FILE *f) = 0;
  virtual bool save(FILE *f) const = 0;
  virtual void save_extra(void) const { };

  int width;
  int height;
};

}  // namespace bhr

#endif
