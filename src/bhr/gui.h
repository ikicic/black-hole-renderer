#ifndef GUI_H
#define GUI_H

namespace bhr {

typedef CartesianVector4<double> _Vector4d;
void render_geodesics(const std::vector<std::pair<_Vector4d, _Vector4d>> &path);

}  // namespace bhr

#endif
