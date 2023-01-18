#ifndef TGA_H
#define TGA_H

#include <bhr/base.h>
#include <bhr/texture.h>

namespace bhr {

bool load_TGA(Image *imaeg, FILE *f);
bool save_TGA(RGBA *image, int width, int height, FILE *f);
bool save_TGA(RGBA *image, int width, int height, const char *filename);
bool save_CSV(RGBd *image, int width, int height, FILE *f);
bool save_CSV(RGBd *image, int width, int height, const char *filename);

}  // namespace bhr

#endif
