#ifndef MOD_H
#define MOD_H

#include <limits>
#include <type_traits>

// http://stackoverflow.com/questions/4633177/c-how-to-wrap-a-float-to-the-interval-pi-pi
template<typename T, typename T2>
T Mod(T x, T2 y) {
  static_assert(!std::numeric_limits<T>::is_exact , "Mod: floating-point type expected");

  if (0. == y)
    return x;

  double m= x - y * floor(x/y);

  if (y > 0)
  {
    if (m>=y)           // Mod(-1e-16             , 360.    ): m= 360.
      return 0;

    if (m<0) {
      if (y+m == y)
        return 0  ; // just in case...
      else
        return y+m; // Mod(106.81415022205296 , _TWO_PI ): m= -1.421e-14
    }
  }
  else                    // modulo range: (y..0]
  {
    if (m<=y)           // Mod(1e-16              , -360.   ): m= -360.
      return 0;

    if (m>0) {
      if (y+m == y)
        return 0  ; // just in case...
      else
        return y+m; // Mod(-106.81415022205296, -_TWO_PI): m= 1.421e-14
    }
  }

  return m;
}

#endif
