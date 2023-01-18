#ifndef TEXTURE_H
#define TEXTURE_H

#include <bhr/base.h>

namespace bhr {

struct Image {
  int width;
  int height;
  int bpp;
  void *data;

  Image() : width(0), height(0), bpp(0), data(NULL) {}
  Image(const char *filename) : Image() {
    load(filename);
  }
  ~Image() {
    free();
  }

  bool load(const char *filename);
  void free(void);

  /* x and y between 0 and width - 1, 0 and height - 1 */
  RGBA get_pixel(int x, int y) const;

  /* x and y are between 0 and 1 */
  RGBA get_pixel_rel(double x, double y) const;
};

Image *load_image_and_cache(const std::string &filename);

struct GUITexture {
  unsigned ID;

  GUITexture() : ID(0) {}
  GUITexture(const Image &image) : ID(0) {
    load(image);
  }

  void bind(void) const;
  void load(const Image &image);
};


void load_image(const char *filename, Image &);
void image_to_texture(const Image &);



class SingleColorTexture {
  RGBd color;
 public:
  struct TexCoord { /* empty */ };
  SingleColorTexture(const RGBd &_color) : color(_color) {}

  inline RGBd get_color(const TexCoord &, const TexCoord &, const TexCoord &,
      const TexCoord &) const {
    return color;
  }
};

class ColorMixTexture {
 public:
  using TexCoord = RGBd;

  inline RGBd get_color(const TexCoord &A, const TexCoord &B, const TexCoord &C,
      const TexCoord &D) const {
    return 0.25 * (A + B + C + D);
  }
};

class DummyTexture {
  struct TexCoord { /* empty */ };
  inline RGBd get_color(const TexCoord &, const TexCoord &, const TexCoord &,
      const TexCoord &) const {
    return RGBd{.5, .5, .5};
  }
};

}  // namespace bhr

#endif
