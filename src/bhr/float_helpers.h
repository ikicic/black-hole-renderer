#ifndef FLOAT_HELPER_H
#define FLOAT_HELPER_H

#include <limits>

namespace bhr {

// Based on https://gist.github.com/alexshtf/eb5128b3e3e143187794
constexpr double constexpr_sqrt(double x) {
    if (!(0.0 <= x && x < std::numeric_limits<double>::infinity()))
        return std::numeric_limits<double>::quiet_NaN();

    double prev = 0.0;
    double curr = x;
    while (prev != curr) {
        prev = curr;
        curr = 0.5 * (curr + x / curr);
    }
    return curr;
}

constexpr float constexpr_sqrt(float x) {
  return static_cast<float>(constexpr_sqrt(static_cast<double>(x)));
}

}  // namespace bhr

#endif
