#ifndef TGA_H
#define TGA_H

#include "texture.h"

bool load_TGA(Image *imaeg, FILE *f);
bool save_TGA(RGBA *image, int width, int height, FILE *f);
bool save_TGA(RGBA *image, int width, int height, const char *filename);
bool save_CSV(RGBd *image, int width, int height, FILE *f);
bool save_CSV(RGBd *image, int width, int height, const char *filename);

#endif
