#ifndef COLORS_H
#define COLORS_H

namespace bhr {

RGBd linear_color_scale(colreal_t low, colreal_t mid, colreal_t high);
RGBd discrete_linear_color_scale(colreal_t low, colreal_t mid, colreal_t high);

}  // namespace bhr

#endif
