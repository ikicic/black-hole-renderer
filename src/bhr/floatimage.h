#ifndef FLOAT_IMAGE_H
#define FLOAT_IMAGE_H

namespace bhr {

RGBd *load_floating_point_image(FILE *f, int width, int height);
bool save_floating_point_image(FILE *f, int width, int height, RGBd *rgb);

}  // namespace bhr

#endif
