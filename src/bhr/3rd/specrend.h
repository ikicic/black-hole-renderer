#ifndef SPECREND_H
#define SPECREND_H

#include "../spectrum.h"

namespace bhr {

/* A colour system is defined by the CIE x and y coordinates of
   its three primary illuminants and the x and y coordinates of
   the white point. */

// Source: http://www.fourmilab.ch/documents/specrend/specrend.c
struct ColourSystem {
  const char *name;             /* Colour system name */
  double xRed, yRed;            /* Red x, y */
  double xGreen, yGreen;        /* Green x, y */
  double xBlue, yBlue;          /* Blue x, y */
  double xWhite, yWhite;        /* White point x, y */
  double gamma;                 /* Gamma correction for system */
};

extern ColourSystem NTSCsystem;
extern ColourSystem EBUsystem;
extern ColourSystem SMPTEsystem;
extern ColourSystem HDTVsystem;
extern ColourSystem CIEsystem;
extern ColourSystem Rec709system;

extern double cie_colour_match[81][3];

RGBd xyz_to_rgb(const ColourSystem &cs, const XYZd &xyz);
bool constrain_rgb(RGBd *rgb);

// Source: http://www.fourmilab.ch/documents/specrend/specrend.c
// Modified to use a functor instead of a function pointer.
template <typename Func>
XYZd spectrum_to_xyz(const Func &spec_intens) {
  /* CIE colour matching functions xBar, yBar, and zBar for
     wavelengths from 380 through 780 nanometers, every 5
     nanometers.  For a wavelength lambda in this range:

        cie_colour_match[(lambda - 380) / 5][0] = xBar
        cie_colour_match[(lambda - 380) / 5][1] = yBar
        cie_colour_match[(lambda - 380) / 5][2] = zBar

  To save memory, this table can be declared as floats
  rather than doubles; (IEEE) float has enough
  significant bits to represent the values. It's declared
  as a double here to avoid warnings about "conversion
  between floating-point types" from certain persnickety
  compilers. */

  int i;
  colreal_t lambda, X = 0, Y = 0, Z = 0;
  for (i = 0, lambda = 380; lambda < colreal_t(780.1); i++, lambda += 5) {
    colreal_t Me = spec_intens(lambda);
    X += Me * cie_colour_match[i][0];
    Y += Me * cie_colour_match[i][1];
    Z += Me * cie_colour_match[i][2];
  }

  // double XYZ = (X + Y + Z);
  // return XYZd{X / XYZ, Y / XYZ, Z / XYZ};
  return XYZd{{X, Y, Z}};
}

}  // namespace bhr

#endif
