#ifndef LINE_H
#define LINE_H

template <typename T, typename Color>
void render_line(T x1, T y1, T x2, T y2,
    Color color, int width, int height, Color *output) {
  // TODO: Better line rendering algorithm.
  auto len = std::sqrt(sqr(x2 - x1) + sqr(y2 - y1));
  for (int i = 0; i <= len; ++i) {
    int x = int(x1 + (x2 - x1) * i / len);
    int y = int(y1 + (y2 - y1) * i / len);
    if (x >= 0 && x < width && y >= 0 && y < height)
      output[y * width + x] += color;
  }
}

template <typename T, typename Color>
void render_arrow(T x, T y, T dx, T dy, T size, T angle,
    Color color, int width, int height, Color *output) {
  T length = std::sqrt(sqr(dx) + sqr(dy));
  T x2 = x + dx;
  T y2 = y + dy;
  if (size > length / 2)
    size = length / 2;
  T ssin = (size / length) * std::sin(angle);
  T scos = (size / length) * std::cos(angle);
  render_line(x, y, x2, y2, color, width, height, output);
  render_line(x2, y2, x2 - dx * scos - dy * ssin, y2 - dy * scos + dx * ssin,
              color, width, height, output);
  render_line(x2, y2, x2 - dx * scos + dy * ssin, y2 - dy * scos - dx * ssin,
              color, width, height, output);
}

#endif
