#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>

#include "base.h"
#include "texture.h"
#include "tga.h"

#define TEXTURE_FILTER_NONE       1
#define TEXTURE_FILTER_BILINEAR   2

#define TEXTURE_FILTER            TEXTURE_FILTER_BILINEAR


bool Image::load(const char *filename) {
  this->free();

  fprintf(stderr, "Loading texture: %s\n", filename);
  FILE *f = fopen(filename, "rb");
  if (f == nullptr) {
    fprintf(stderr, "Loading failed!\n");
    return false;
  }

  bool result = load_TGA(this, f);
  fclose(f);
  if (!result) fprintf(stderr, "Error loading texture %s\n", filename);
  return result;
}

void Image::free(void) {
  ::free(data);
  data = NULL;
}


RGBA Image::get_pixel(int x, int y) const {
  unsigned char *p = (unsigned char *)data + (y * width + x) * (bpp / 8);
  return RGBA(p[0], p[1], p[2]);
}

inline Vector<double, 3> _get_pixel(const Image *image, int x, int y) {
  unsigned char *p = (unsigned char *)image->data + (y * image->width + x) * (image->bpp / 8);
  return {{double(p[0]), double(p[1]), double(p[2])}};
}

RGBA Image::get_pixel_rel(double x, double y) const {
#if TEXTURE_FILTER == TEXTURE_FILTER_NONE
  return get_pixel(int(x * (width - 1)), int(y * (height - 1)));
#elif TEXTURE_FILTER == TEXTURE_FILTER_BILINEAR
  x *= width - 1;
  y *= height - 1;
  int X = int(x); x -= X;
  int Y = int(y); y -= Y;
  if (X == width - 1 || Y == height - 1) {
    fprintf(stderr, "sranje %d %d  %lf %lf\n", X, Y, x, y);
    exit(0);
  }
  auto result =
    ((1 - x) * (1 - y)) * _get_pixel(this, X, Y) +
    ((1 - x) * y) * _get_pixel(this, X, Y + 1) +
    (x * (1 - y)) * _get_pixel(this, X + 1, Y) +
    (x * y) * _get_pixel(this, X + 1, Y + 1);

  return RGBA(
      (unsigned char)(int)result[0],
      (unsigned char)(int)result[1],
      (unsigned char)(int)result[2]
    );
#endif
}

Image *load_image_and_cache(const std::string &filename) {
  static std::map<std::string, Image *> cache;
  auto it = cache.find(filename);
  if (it != cache.end())
    return it->second;

  auto iter = cache.emplace(filename, new Image(filename.c_str())).first;
  if (iter->second->data == nullptr) {
    delete iter->second;
    iter->second = nullptr;
  }

  return iter->second;
}
