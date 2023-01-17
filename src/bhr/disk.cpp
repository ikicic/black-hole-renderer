#include <bhr/base.h>

#if RENDER_DISK
#include <bhr/disk.h>

ShakuraSunyaevDisk *_shakura_sunyaev = nullptr;
double shakura_sunyaev_height(double r) {
  return _shakura_sunyaev->get_height(r);
}

#if DISK_RELIEF_TEXTURE
Image disk_relief_tex("images/relief.tga");

bool load_disk_relief_texture(void) {
  if (disk_relief_tex.data != nullptr)
    return true;
  return disk_relief_tex.load("images/relief.tga");
}
#endif

#endif
