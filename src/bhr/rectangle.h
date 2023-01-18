#ifndef RECTANGLE_H
#define RECTANGLE_H

template <typename T>
inline T ccw(T x1, T y1, T x2, T y2, T x3, T y3) {
  return x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2);
}

template <typename T, typename TexCoord, typename Color, typename Texture>
void render_textured_axes_aligned_rectangle(
    T x1, T y1,
    T x2, T y2,
    TexCoord t1,
    TexCoord t2,
    TexCoord t3,
    TexCoord t4,
    const Texture &texture,
    int width, int height,
    Color *output) {

  /*  (x1, y1).................
   *     .                    .
   *     .                    .
   *     ................. (x2, y2)
   *
   * (tx1, ty1)...........(tx2, ty2)
   *     .                    .
   *     .                    .
   *     .                    .
   * (tx4, ty4)...........(tx3, ty3)
   */

  if (x2 < x1) {
    std::swap(x1, x2);
    std::swap(t1, t2);
    std::swap(t3, t4);
  }
  if (y2 < y1) {
    std::swap(y1, y2);
    std::swap(t1, t4);
    std::swap(t2, t3);
  }

  if (x2 <= 0 || x1 >= width) return;
  if (y2 <= 0 || y1 >= height) return;

  if (x1 < 0) {
    double c = (-x1 / (x2 - x1));
    t1 = (1 - c) * t1 + c * t2;
    t4 = (1 - c) * t4 + c * t3;
    x1 = 0;
  }
  if (x2 > width) {
    double c = (x2 - width) / (x2 - x1);
    t2 = (1 - c) * t2 + c * t1;
    t3 = (1 - c) * t3 + c * t4;
    x2 = width;
  }
  if (y1 < 0) {
    double c = (-y1 / (y2 - y1));
    t1 = (1 - c) * t1 + c * t4;
    t2 = (1 - c) * t2 + c * t3;
    y1 = 0;
  }
  if (y2 > height) {
    double c = (y2 - height) / (y2 - y1);
    t4 = (1 - c) * t4 + c * t1;
    t3 = (1 - c) * t3 + c * t2;
    y2 = height;
  }

  int X1 = (int)x1;
  int Y1 = (int)y1;
  int X2 = (int)x2;
  int Y2 = (int)y2;

  if (X1 >= width) { X1 = X2 = width - 1; } else if (X2 >= width) X2 = width - 1;
  if (Y1 >= height) { Y1 = Y2 = height - 1; } else if (Y2 >= height) Y2 = height - 1;

  // http://www.uio.no/studier/emner/matnat/ifi/INF4360/h09/undervisningsmateriale/bc_lecture.pdf
  // http://math.arizona.edu/~agillette/research/ccomOct11.pdf
  auto wachspress_mix = [t1, t2, t3, t4, x1, y1, x2, y2](T x, T y) {
    T A1 = ccw(x1, y1, x2, y1, x, y);
    T A2 = ccw(x2, y1, x2, y2, x, y);
    T A3 = ccw(x2, y2, x1, y2, x, y);
    T A4 = ccw(x1, y2, x1, y1, x, y);

    T w1 = A2 * A3;
    T w2 = A3 * A4;
    T w3 = A4 * A1;
    T w4 = A1 * A2;
    T invw = T(1) / (w1 + w2 + w3 + w4);

    return (w1 * invw) * t1
         + (w2 * invw) * t2
         + (w3 * invw) * t3
         + (w4 * invw) * t4;
  };

  auto ADD = [output, &wachspress_mix, &texture, &width](
      int X, int Y, T xa, T ya, T xb, T yb) {

    // Wachspress Coordinates.
    TexCoord tf1 = wachspress_mix(xa, ya);
    TexCoord tf2 = wachspress_mix(xb, ya);
    TexCoord tf3 = wachspress_mix(xb, yb);
    TexCoord tf4 = wachspress_mix(xa, yb);

    Color color = texture.get_color(tf1, tf2, tf3, tf4);
    output[Y * width + X] += ((xb - xa) * (yb - ya)) * color;
  };

  // TODO: add x1, y1 to all T coords.
  if (X1 == X2) {
    if (Y1 == Y2) {
      ADD(X1, Y1, x1, y1, x2, y2);
    } else {
      ADD(X1, Y1, x1, y1, x2, Y1 + 1);
      for (int i = Y1 + 1; i < Y2; ++i)
        ADD(X1, i, x1, i, x2, i + 1);
      ADD(X1, Y2, x1, Y2, x2, y2);
    }
  } else {
    if (Y1 == Y2) {
      ADD(X1, Y1, x1, y1, X1 + 1, y2);
      for (int j = X1 + 1; j < X2; ++j)
        ADD(j, Y1, j, y1, j + 1, y2);
      ADD(X2, Y1, X2, y1, x2, y2);
    } else {
      ADD(X1, Y1, x1, y1, X1 + 1, Y1 + 1);
      for (int j = X1 + 1; j < X2; ++j)
        ADD(j, Y1, j, y1, j + 1, Y1 + 1);
      ADD(X2, Y1, X2, y1, x2, Y1 + 1);
      for (int i = Y1 + 1; i < Y2; ++i) {
        ADD(X1, i, x1, i, X1 + 1, i + 1);
        for (int j = X1 + 1; j < X2; ++j)
          ADD(j, i, j, i, j + 1, i + 1);
        ADD(X2, i, X2, i, x2, i + 1);
      }
      ADD(X1, Y2, x1, Y2, X1 + 1, y2);
      for (int j = X1 + 1; j < X2; ++j)
        ADD(j, Y2, j, Y2, j + 1, y2);
      ADD(X2, Y2, X2, Y2, x2, y2);
    }
  }
#undef ADD
};

template <typename T, typename Color>
void render_axes_aligned_rectangle(T x1, T y1, T x2, T y2,
    Color color, int width, int height, Color *output) {
  /* A pixel (j, i) is represented by x in [j, j + 1> and y in [i, i + 1>. */

  if (x2 < x1) std::swap(x1, x2);
  if (y2 < y1) std::swap(y1, y2);

  x1 = clamp(x1, T(0), (T)width);
  y1 = clamp(y1, T(0), (T)height);
  x2 = clamp(x2, T(0), (T)width);
  y2 = clamp(y2, T(0), (T)height);

  int X1 = (int)x1;
  int Y1 = (int)y1;
  int X2 = (int)x2;
  int Y2 = (int)y2;

  if (X1 >= width) { X1 = X2 = width - 1; } else if (X2 >= width) X2 = width - 1;
  if (Y1 >= height) { Y1 = Y2 = height - 1; } else if (Y2 >= height) Y2 = height - 1;

#define ADD(X, Y, factor) output[(Y) * width + (X)] += (factor) * color
  if (X1 == X2) {
    if (Y1 == Y2) {
      ADD(X1, Y1, (y2 - y1) * (x2 - x1));
    } else {
      ADD(X1, Y1, (Y1 + 1 - y1) * (x2 - x1));
      for (int i = Y1 + 1; i < Y2; ++i)
        ADD(X1, i, x2 - x1);
      ADD(X1, Y2, (y2 - Y2) * (x2 - x1));
    }
  } else {
    if (Y1 == Y2) {
      ADD(X1, Y1, (X1 + 1 - x1) * (y2 - y1));
      for (int j = X1 + 1; j < X2; ++j)
        ADD(j, Y1, y2 - y1);
      ADD(X2, Y1, (x2 - X2) * (y2 - y1));
    } else {
      ADD(X1, Y1, (X1 + 1 - x1) * (Y1 + 1 - y1));
      for (int j = X1 + 1; j < X2; ++j)
        ADD(j, Y1, Y1 + 1 - y1);
      ADD(X2, Y1, (x2 - X2) * (Y1 + 1 - y1));
      for (int i = Y1 + 1; i < Y2; ++i) {
        ADD(X1, i, X1 + 1 - x1);
        for (int j = X1 + 1; j < X2; ++j)
          ADD(j, i, 1);
        ADD(X2, i, x2 - X2);
      }
      ADD(X1, Y2, (X1 + 1 - x1) * (y2 - Y2));
      for (int j = X1 + 1; j < X2; ++j)
        ADD(j, Y2, y2 - Y2);
      ADD(X2, Y2, (x2 - X2) * (y2 - Y2));
    }
  }
#undef ADD
}

template <typename T, typename Color>
void render_axes_aligned_rectangle_border(T x1, T y1, T x2, T y2,
    Color color,
    double inner, double outer,
    int width, int height, Color *output) {
  if (x2 < x1) std::swap(x1, x2);
  if (y2 < y1) std::swap(y1, y2);

  // Top, Bottom, Left, Right.
#define RENDER(_x1, _y1, _x2, _y2) render_axes_aligned_rectangle( \
      _x1, _y1, _x2, _y2, color, width, height, output);
  // RENDER(x1 - outer, y1 - outer, x2 + outer, y1 + inner);
  // RENDER(x1 - outer, y2 + outer, x2 + outer, y2 - inner);
  // RENDER(x1 - outer, y1 + inner, x1 + inner, y2 - inner);
  // RENDER(x2 - outer, y1 + inner, x2 + inner, y2 - inner);
  RENDER(x1 - outer, y1 - outer, x2 - inner, y1 + inner);
  RENDER(x1 - outer, y1 + inner, x1 + inner, y2 - inner);
#undef RENDER
}

#endif
