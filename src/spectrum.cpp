#include "base.h"
#include "spectrum.h"
#include "utility.h"
#include "3rd/specrend.h"

#include <cassert>
#include <cmath>
#include <vector>

#define MIN_TEMPERATURE   1.
#define MAX_TEMPERATURE   100000.
#define TEMPERATURE_POINTS  1000

std::vector<RGBd> bb_rgb_colors;

class BlackBodySpectrum {
  /* Black body spectrum:
   *
   *                2 h c^2                1
   * I(lambda, T) = --------  ---------------------------
   *                lambda^5  exp(h c / (lambda k T)) - 1
   */
 private:
  static constexpr double h = 6.626070040e-34;
  static constexpr double c = 299792458.;
  static constexpr double k = 1.38064852e-23;

  double hc_over_kT;
 public:

  BlackBodySpectrum(double T) : hc_over_kT((h * c / k) / T) {}

  inline double operator()(double wavelength_nm) const {
    const double wavelength = wavelength_nm * 1e-9;
    return (2 * h * c * c)
        / (std::pow(wavelength, 5) * std::expm1(hc_over_kT / wavelength));
  }
};

void calculate_spectrum_colors(void) {
  bb_rgb_colors.reserve(TEMPERATURE_POINTS);
  for (int i = 0; i < TEMPERATURE_POINTS; ++i) {
    double T = MIN_TEMPERATURE
        + (MAX_TEMPERATURE - MIN_TEMPERATURE) * i / (TEMPERATURE_POINTS - 1);
    XYZd xyz = spectrum_to_xyz(BlackBodySpectrum(T));
    bb_rgb_colors.push_back(xyz_to_rgb(HDTVsystem, xyz));
  }
}

RGBd get_black_body_color(double T) {
  assert(!bb_rgb_colors.empty());
  if (T < MIN_TEMPERATURE)
    return bb_rgb_colors[0];
  double k = (TEMPERATURE_POINTS - 1)
      * ((T - MIN_TEMPERATURE) / (MAX_TEMPERATURE - MIN_TEMPERATURE));
  int K = (int)k;
  if (K >= TEMPERATURE_POINTS - 1)
    return bb_rgb_colors.back();

  double high = k - K;
  return (1 - high) * bb_rgb_colors[K] + high * bb_rgb_colors[k + 1];
}

RGBd get_black_body_color(double local_temperature, double freq_factor) {
  /* observer_frequency == freq_factor * local_frequency.  */

  // Yes, no other factors.
  return get_black_body_color(local_temperature * freq_factor);
}

RGBA to_RGBAub(const RGBd &rgb) {
  return RGBA(
      clamp((int)(255 * rgb[0]), 0, 255),
      clamp((int)(255 * rgb[1]), 0, 255),
      clamp((int)(255 * rgb[2]), 0, 255)
  );
}
