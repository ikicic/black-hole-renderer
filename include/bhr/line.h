#ifndef LINE_H
#define LINE_H

template <typename _T, typename _Color>
void render_line(_T x1, _T y1, _T x2, _T y2,
    _Color color, int width, int height, _Color *output) {
  // TODO: Better line rendering algorithm.
  auto len = std::sqrt(sqr(x2 - x1) + sqr(y2 - y1));
  for (int i = 0; i <= len; ++i) {
    int x = int(x1 + (x2 - x1) * i / len);
    int y = int(y1 + (y2 - y1) * i / len);
    if (x >= 0 && x < width && y >= 0 && y < height)
      output[y * width + x] += color;
  }
}

template <typename _T, typename _Color>
void render_arrow(_T x, _T y, _T dx, _T dy, _T size, _T angle,
    _Color color, int width, int height, _Color *output) {
  _T length = std::sqrt(sqr(dx) + sqr(dy));
  _T x2 = x + dx;
  _T y2 = y + dy;
  if (size > length / 2)
    size = length / 2;
  _T ssin = (size / length) * std::sin(angle);
  _T scos = (size / length) * std::cos(angle);
  render_line(x, y, x2, y2, color, width, height, output);
  render_line(x2, y2, x2 - dx * scos - dy * ssin, y2 - dy * scos + dx * ssin,
              color, width, height, output);
  render_line(x2, y2, x2 - dx * scos + dy * ssin, y2 - dy * scos - dx * ssin,
              color, width, height, output);
}

#endif
