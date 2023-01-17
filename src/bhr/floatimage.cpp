#include <bhr/base.h>
#include <bhr/floatimage.h>


RGBd *load_floating_point_image(FILE *f, int width, int height) {
  unsigned int size = width * height;
  RGBd *rgbd = new RGBd[size];
  RGBf *rgbf = new RGBf[size];
  if (rgbd == nullptr || rgbf == nullptr) {
    fprintf(stderr, "Couldn't allocate %dx%d RGBd and RGBf (%uB).\n",
        width, height, (unsigned int)(size * sizeof(RGBd) + sizeof(RGBf)));
    return nullptr;
  }
  if (fread(rgbf, sizeof(RGBf), size, f) != size) {
    delete []rgbd;
    delete []rgbf;
    return nullptr;
  }
  for (unsigned int k = 0; k < size; ++k)
    rgbd[k] = RGBd(rgbf[k]);
  delete []rgbf;
  return rgbd;
}

bool save_floating_point_image(FILE *f, int width, int height, RGBd *rgbd) {
  unsigned int size = width * height;
  RGBf *rgbf = new RGBf[size];
  if (rgbf == nullptr) {
    fprintf(stderr, "Couldn't allocate %uB for double->float conversion.\n",
        (unsigned int)(size * sizeof(RGBf)));
    return false;
  }
  for (unsigned int k = 0; k < size; ++k)
    rgbf[k] = RGBf(rgbd[k]);
  unsigned int saved = fwrite(rgbf, sizeof(RGBf), size, f);
  delete []rgbf;
  return saved == size;
}
