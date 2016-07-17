#ifndef FLOAT_IMAGE_H
#define FLOAT_IMAGE_H

RGBd *load_floating_point_image(FILE *f, int width, int height);
bool save_floating_point_image(FILE *f, int width, int height, RGBd *rgb);

#endif
