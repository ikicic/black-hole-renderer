#ifndef SETTINGS_H
#define SETTINGS_H

#include <bhr/coordinate.h>

namespace bhr {

enum SpacetimeEnum {
  SPACETIME_FLAT = 1,
  SPACETIME_SCHWARZSCHILD = 2,
  SPACETIME_KERR = 3,
  SPACETIME_REISSNER_NORDSTROM = 4,
};

class Settings {
 public:
  // Camera
  Vector3 camera_position_cart;
  Vector3 camera_to_cart;
  Vector3 camera_up_cart;
  real_t fovy;
  real_t ver_range_km;
  bool ortho;

  // preprocessed data
  std::string input_filename;

  // output image
  int width;
  int height;
  std::string output_geodesics;
  std::string output_float;
  std::string output_image;
  bool cache;
  bool horizontal_flip;

  // spacetime
  SpacetimeEnum type;
  real_t M;  // Mass.
  real_t a;  // J / M.
  real_t Q;  // Charge.

  // field
  Vector3 field;    // B (in teslas), theta and phi (in degrees)

  // debug
  int debug;

  // rendering
  std::string sky_image;
  real_t dlambda;
  int threads;  // Not included in the hash.
  bool recursive;
  int max_extra_recursive_depth;  // max =~ log2(max(width, height)) + extra.

  bool check(void) const;
  void set_default(void);
  bool read_cli(int agrc, char **argv);

  std::string get_geodesics_hash(void) const;
  std::string get_float_hash(void) const;
  std::string get_image_hash(void) const;
  std::string get_auto_geodesics_filename(void) const;
  std::string get_float_filename(void) const;
  std::string get_image_filename(void) const;
};

void check_settings(const Settings &settings);
void default_settings(Settings &settings);
bool cli_to_settings(int argc, char **argv, Settings &settings);

std::string get_settings_geodesics_hash(const Settings &settings);
std::string get_settings_image_hash(const Settings &settings);
std::string generate_geodesics_filename(const Settings &settings);
std::string generate_image_filename(const Settings &settings);

}  // namespace bhr

#endif
