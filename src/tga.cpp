#include "base.h"
#include "tga.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

bool load_TGA(Image *image, FILE *f) {
  char file_header[12];
  const char tga_header[12] = {0,0,2,0,0,0,0,0,0,0,0,0};
  unsigned char header[6];

  fread(file_header, 12, 1, f);
  if (memcmp(file_header, tga_header, 12)) {
    fprintf(stderr, "Unsupported TGA header\n");
    return false;
  }

  fread(header, 6, 1, f);

  image->width = ((unsigned)header[1] << 8) | (unsigned)header[0];
  image->height = ((unsigned)header[3] << 8) | (unsigned)header[2];
  fprintf(stderr, "Texture size: %d %d\n", image->width, image->height);
  image->bpp = header[4];

  int size = (image->width * image->height * image->bpp + 7) / 8;
  image->data = malloc(size);
  if (image->data == NULL)
    return false;

  int result = fread(image->data, size, 1, f);
  if (result != 1 || ferror(f)) {
    fprintf(stderr, "Error reading %d bytes!\n", size);
    free(image->data);
    image->data = nullptr;
    return false;
  }

  // SWAP RB
  return true;
}

bool save_TGA(RGBA *image, int width, int height, FILE *f) {
  unsigned char header[18] = {0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  header[12] = width         & 0xFF;
  header[13] = (width >> 8)  & 0xFF;
  header[14] = height        & 0xFF;
  header[15] = (height >> 8) & 0xFF;
  header[16] = 32;
  header[17] = 32;  // Vertical flip.

  if (fwrite(header, 18, 1, f) != 1) return false;

  const int size = width * height;
  RGBA *buffer = new RGBA[size];
  if (buffer == nullptr) {
    fprintf(stderr, "Allocation of %db failed!\n", size);
    return false;
  }
  // Convert RGBA to BGRA.
  for (int k = 0; k < size; ++k)
    buffer[k] = image[k].get_BGRA();
  if (fwrite(buffer, size * 4, 1, f) != 1) {
    delete []buffer;
    return false;
  }
  delete []buffer;
  return true;
}

bool save_TGA(RGBA *image, int width, int height, const char *filename) {
  FILE *f = fopen(filename, "wb");
  if (f == nullptr) {
    fprintf(stderr, "save_TGA: File %s not found!\n", filename);
    return false;
  }

  bool result = save_TGA(image, width, height, f);
  fclose(f);
  return result;
}

bool save_CSV(RGBd *image, int width, int height, FILE *f) {
  for (int i = 0; i < height; ++i)
    for (int j = 0; j < width; ++j) {
      int result = fprintf(f, j == width - 1 ? "%lg,%lg,%lg\n" : "%lg,%lg,%lg,",
          image[i * width + j][0],
          image[i * width + j][1],
          image[i * width + j][2]
      );
      if (result < 0)
        return false;
    }
  return true;
}

bool save_CSV(RGBd *image, int width, int height, const char *filename) {
  FILE *f = fopen(filename, "w");
  if (f == nullptr) {
    fprintf(stderr, "save_CSV: Error opening file %s!\n", filename);
    return false;
  }

  bool result = save_CSV(image, width, height, f);
  fclose(f);
  return result;

}
