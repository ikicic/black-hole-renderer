#include "base.h"
#include "settings.h"
#include "parameters.h"

#include "3rd/sha1.h"

#include <sstream>
#include <thread>

bool Settings::check(void) const {
  if (type == SPACETIME_KERR && a > M) {
    fprintf(stderr, "Warning: a > M! (a=%lf, M=%lf)\n",
        (double)a, (double)M);
  }
  return true;
}

void Settings::set_default(void) {
  constexpr double r_s = _BLACK_HOLE_r_S;
#if RENDER_DISK == DISK_KERTAP
  constexpr double theta = 75 * M_PI / 180;
  constexpr double phi = 180 * M_PI / 180;
#elif RENDER_DISK == DISK_DUMMY
  constexpr double theta = 80 * M_PI / 180;
  constexpr double phi = 180 * M_PI / 180;
#elif RENDER_DISK == DISK_SHAKURA
  constexpr double theta = 75 * M_PI / 180;
  constexpr double phi = 0 * M_PI / 180;
#else
  constexpr double theta = 90 * M_PI / 180;
  constexpr double phi = 180 * M_PI / 180;
#endif
  // constexpr double phi = 0 * M_PI / 180;
  constexpr double dist = 0.9 * MAX_r;
  camera_position_cart = {{
    dist * std::sin(theta) * std::cos(phi),
    dist * std::sin(theta) * std::sin(phi),
    dist * std::cos(theta)
  }};
  // camera_position_cart = {{0, 100 * r_s, 0}};
  // camera_to_cart = {0.2, 0, 0};
  // camera_to_cart = {{0.00001, 0.00009, 0.0000013}};
  camera_to_cart = {{0, 0, 0}};
  camera_up_cart = {{0, 0, 1}};
  fovy = 0.09;
  ver_range_km = r_s * 400 / UNIT_km;
  ortho = false;

  input_filename = "";

  width = 200;
  height = 200;
  output_float = "";
  output_image = "auto";
  cache = true;
  horizontal_flip = false;

  // type = SPACETIME_FLAT;
  type = SPACETIME_KERR;
  // type = SPACETIME_SCHWARZSCHILD;
  // type = SPACETIME_REISSNER_NORDSTROM;
  M = BLACK_HOLE_M;
  a = BLACK_HOLE_a;
  Q = BLACK_HOLE_Q;

  field = Vector3{2e11, 0., 0.};

  debug = 0;

  sky_image = "images/sky_high.tga";
  dlambda = 0.100001;
  threads = std::max((int)std::thread::hardware_concurrency(), 1);
  recursive = false;
  max_extra_recursive_depth = 0;
}

static void _read_value(int &value, const char *input) {
  sscanf(input, "%d", &value);
}
static void _read_value(float &value, const char *input) {
  sscanf(input, "%f", &value);
}
static void _read_value(double &value, const char *input) {
  sscanf(input, "%lf", &value);
}
static void _read_value(long double &value, const char *input) {
  sscanf(input, "%Lf", &value);
}
static void _read_value(Vector<float, 3> &value, const char *input) {
  sscanf(input, "%f,%f,%f", &value[0], &value[1], &value[2]);
}
static void _read_value(Vector<double, 3> &value, const char *input) {
  sscanf(input, "%lf,%lf,%lf", &value[0], &value[1], &value[2]);
}
static void _read_value(Vector<long double, 3> &value, const char *input) {
  sscanf(input, "%Lf,%Lf,%Lf", &value[0], &value[1], &value[2]);
}
static void _read_value(std::string &value, const char *input) {
  value = input;
}

template<typename _T> static int _read(
    int argc, char **argv, int &i,
    const std::string &cur, const std::string &name, _T &value) {
  if (cur == "--" + name) {
    if (i == argc - 1) {
      fprintf(stderr, "Missing value for %s\n", cur.c_str());
      return -1;
    } else {
      _read_value(value, argv[i + 1]);
      i += 1;
      return 1;
    }
  }
  return 0;
}

static int _read(
    int /*argc*/, char ** /*argv*/, int &/*i*/,
    const std::string &cur, const std::string &name, bool &value) {
  if (cur == "--" + name)
    return value = true, 1;
  if (cur == "--no-" + name)
    return value = false, 1;
  return 0;
}

static int _read(
    int /*argc*/, char ** /*argv*/, int &/*i*/,
    const std::string &cur, const std::string &/*name*/,
    SpacetimeEnum &value) {
  if (cur == "flat") value = SPACETIME_FLAT;
  else if (cur == "schwarzschild") value = SPACETIME_SCHWARZSCHILD;
  else if (cur == "kerr") value = SPACETIME_KERR;
  else return 0;
  return 1;
}

#define READ(name) { \
    int _dummy = _read(argc, argv, i, cur, std::string(#name), name); \
    if (_dummy == 1) continue; \
    if (_dummy == -1) return false; \
  }
#define REJECT(name) { \
    if (std::string(argv[i]) == #name) { \
      fprintf(stderr, "Argument \"%s\" is disabled.\n", #name); \
      return false; \
    } \
  }
bool Settings::read_cli(int argc, char **argv) {
  int _width = width;
  int _height = height;
  int _aspect_x = -1, _aspect_y = 1;
  width = -1;
  height = -1;

  for (int i = 1; i < argc; ++i) {
    std::string cur = argv[i];
    if (cur == "--aspect") {
      if (i < argc - 1) {
        sscanf(argv[++i], "%d:%d", &_aspect_x, &_aspect_y);
        continue;
      } else {
        fprintf(stderr, "Missing value for aspect ratio\n");
        return false;
      }
    }

    READ(camera_position_cart);
    READ(camera_to_cart);
    READ(camera_up_cart);
    READ(fovy);
    READ(ver_range_km);
    READ(ortho);

    READ(input_filename);

    READ(width);
    READ(height);
    READ(output_float);
    READ(output_image);
    READ(cache);
    READ(horizontal_flip);

    READ(type);
#if PREDEFINED_PARAMS
    REJECT(M);
    REJECT(a);
    REJECT(Q);
#else
    READ(M);
    READ(a);
    READ(Q);
#endif

    READ(field);

    READ(debug);

    READ(sky_image);
    READ(dlambda);
    READ(threads);
    READ(recursive);
    READ(max_extra_recursive_depth);

    fprintf(stderr, "Unknown argument: %s\n", argv[i]);
    return false;
  }

  if (_aspect_x != -1 && _aspect_y != -1) {
    if (width != -1) {
      height = width * _aspect_y / _aspect_x;
    } else if (height != -1) {
      width = height * _aspect_x / _aspect_y;
    } else {
      fprintf(stderr, "Resolution not specified\n");
      return false;
    }
  } else {
    if (width == -1) width = _width;
    if (height == -1) height = _height;
  }

  return true;
}


#define DELIM "#"
#define ADD_VECTOR3(v) (s << v[0] << DELIM << v[1] << DELIM << v[2] << DELIM)
std::string Settings::get_geodesics_hash(void) const {
  /* Includes all variables affecting geodesic. */
  std::stringstream s;
  s.precision(30);
  ADD_VECTOR3(camera_position_cart);
  ADD_VECTOR3(camera_to_cart);
  ADD_VECTOR3(camera_up_cart);
  s << ortho << DELIM;
  if (ortho) {
    s << ver_range_km << DELIM;
  } else {
    s << fovy << DELIM;
  }

  s << width << DELIM;
  s << height << DELIM;
  s << horizontal_flip << DELIM;

  s << type << DELIM;
  if (type != SPACETIME_FLAT) {
    s << M / UNIT_kg << DELIM;
  }
  if (type == SPACETIME_KERR) {
    s << a << DELIM;
  }
  if (type == SPACETIME_REISSNER_NORDSTROM) {
    s << Q << DELIM;
  }
  s << field << DELIM;
  s << RENDER_DISK << DELIM;

  s << dlambda << DELIM;
  s << recursive << DELIM;
  s << max_extra_recursive_depth << DELIM;

#if RENDER_DISK
  s << "Disk" << DELIM;
  s << INNER_RADIUS << DELIM;
  s << OUTER_RADIUS << DELIM;
#endif
  return sha1(s.str());
}

std::string Settings::get_float_hash(void) const {
  std::stringstream s;
  if (input_filename.empty())
    s << get_geodesics_hash() << DELIM;

  s << sky_image << DELIM;
  return sha1(s.str());
}

std::string Settings::get_image_hash(void) const {
  return get_float_hash();
}

std::string Settings::get_auto_geodesics_filename(void) const {
  return "output/raw/geo" + get_geodesics_hash() + ".geo";
}

std::string Settings::get_float_filename(void) const {
  return output_float == "auto"
      ? "output/float/float" + get_float_hash() + ".flt"
      : output_float;
}

std::string Settings::get_image_filename(void) const {
  return output_image == "auto"
      ? "output/img" + get_image_hash() + ".tga"
      : output_image;
}

#undef DELIM
#undef ADD_VECTOR3
