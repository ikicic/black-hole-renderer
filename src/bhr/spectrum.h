// TODO: rename to colors.h
#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <bhr/base.h>

namespace bhr {

void calculate_spectrum_colors(void);
RGBd get_black_body_color(double temperature);
RGBd get_black_body_color(double local_temperature, double freq_factor);
RGBA to_RGBAub(const RGBd &rgb);

}  // namespace bhr

#endif
