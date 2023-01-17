#ifndef MATRIX_H
#define MATRIX_H

#include <bhr/utility.h>

template <typename _T, int _N>
struct Matrix {
  _T v[_N][_N];

  inline const _T * operator[](int k) const {
    return v[k];
  }
  inline _T * operator[](int k) {
    return v[k];
  }

  template <template <typename> class _Vector>
  friend inline _Vector<_T> operator*(const Matrix &mat,
                                      const _Vector<_T> &vector) {
    _Vector<_T> result;
    for (int i = 0; i < _N; ++i) {
      _T tmp = _T();
      for (int j = 0; j < _N; ++j)
        tmp += mat.v[i][j] * vector[j];
      result[i] = tmp;
    }
    return result;
  }

  friend inline Matrix operator*(const Matrix &A, const Matrix &B) {
    Matrix result = Matrix();
    for (int k = 0; k < _N; ++k) {
      for (int i = 0; i < _N; ++i) {
        const _T &tmp = A.v[i][k];
        for (int j = 0; j < _N; ++j)
          result.v[i][j] += tmp * B.v[k][j];
      }
    }
    return result;
  }

#define FOR(i, j) for (int i = 0; i < _N; ++i) for (int j = 0; j < _N; ++j)
  friend inline Matrix operator+(const Matrix &A, const Matrix &B) {
    Matrix result;
    FOR(i, j)
      result.v[i][j] = A.v[i][j] + B.v[i][j];
    return result;
  }
  friend inline Matrix operator-(const Matrix &A, const Matrix &B) {
    Matrix result;
    FOR(i, j)
      result.v[i][j] = A.v[i][j] - B.v[i][j];
    return result;
  }
  friend inline Matrix operator*(const _T &a, const Matrix &mat) {
    Matrix result;
    FOR(i, j)
      result.v[i][j] = a * mat.v[i][j];
    return result;
  }

  template <template <typename> class _Vector>
  friend inline _Vector<_T> operator*(const _Vector<_T> &vector,
                                      const Matrix &mat) {
    _Vector<_T> result = _Vector<_T>();
    FOR(i, j)
      result[j] += vector[i] * mat.v[i][j];
    return result;
  }

  inline Matrix &operator+=(const Matrix &mat) {
    FOR(i, j)
      v[i][j] += mat.v[i][j];
    return *this;
  }

  friend inline auto numerical_sqr_distance(const Matrix &A, const Matrix &B) {
    _T result = _T();
    FOR(i, j)
      result += numerical_sqr_distance(A.v[i][j], B.v[i][j]);
    return result;
  }

  friend inline void mult(Matrix *self, const auto &c) {
    FOR(i, j)
      mult(&self->v[i][j], c);
  }

  friend inline void mult_add(Matrix *self, const auto &c, const Matrix &A) {
    FOR(i, j)
      mult_add(&self->v[i][j], c, A.v[i][j]);
  }

  friend inline void set_and_mult_add(
      Matrix *self, const Matrix &A, const auto &c, const Matrix &B) {
    FOR(i, j)
      set_and_mult_add(&self->v[i][j], A.v[i][j], c, B.v[i][j]);
  }
#undef FOR

  friend std::ostream& operator<<(std::ostream& stream, const Matrix &A) {
    for (int i = 0; i < _N; ++i)
      for (int j = 0; j < _N; ++j)
        stream << A.v[i][j] << (j < _N - 1 ? ' ' : '\n');
    return stream << '\n';
  }
};
template <typename _T> using Matrix3 = Matrix<_T, 3>;
template <typename _T> using Matrix4 = Matrix<_T, 4>;


// http://stackoverflow.com/questions/2922690/calculating-an-nxn-matrix-determinant-in-c-sharp/2922905#2922905
template <typename _T>
_T determinant4(const Matrix4<_T> &mat) {
  const _T *m = (const _T *)mat.v;
  return m[12] * m[9]  * m[6]  * m[3]   -  m[8] * m[13] * m[6]  * m[3]   -
         m[12] * m[5]  * m[10] * m[3]   +  m[4] * m[13] * m[10] * m[3]   +
         m[8]  * m[5]  * m[14] * m[3]   -  m[4] * m[9]  * m[14] * m[3]   -
         m[12] * m[9]  * m[2]  * m[7]   +  m[8] * m[13] * m[2]  * m[7]   +
         m[12] * m[1]  * m[10] * m[7]   -  m[0] * m[13] * m[10] * m[7]   -
         m[8]  * m[1]  * m[14] * m[7]   +  m[0] * m[9]  * m[14] * m[7]   +
         m[12] * m[5]  * m[2]  * m[11]  -  m[4] * m[13] * m[2]  * m[11]  -
         m[12] * m[1]  * m[6]  * m[11]  +  m[0] * m[13] * m[6]  * m[11]  +
         m[4]  * m[1]  * m[14] * m[11]  -  m[0] * m[5]  * m[14] * m[11]  -
         m[8]  * m[5]  * m[2]  * m[15]  +  m[4] * m[9]  * m[2]  * m[15]  +
         m[8]  * m[1]  * m[6]  * m[15]  -  m[0] * m[9]  * m[6]  * m[15]  -
         m[4]  * m[1]  * m[10] * m[15]  +  m[0] * m[5]  * m[10] * m[15];
}

/* http://stackoverflow.com/questions/1148309/inverting-a-4x4-matrix */
template <typename _T>
_T gluInvertMatrix(const _T m[16], _T invOut[16]) {
  _T inv[16];

  inv[0] = m[5]  * m[10] * m[15] -
    m[5]  * m[11] * m[14] -
    m[9]  * m[6]  * m[15] +
    m[9]  * m[7]  * m[14] +
    m[13] * m[6]  * m[11] -
    m[13] * m[7]  * m[10];

  inv[4] = -m[4]  * m[10] * m[15] +
    m[4]  * m[11] * m[14] +
    m[8]  * m[6]  * m[15] -
    m[8]  * m[7]  * m[14] -
    m[12] * m[6]  * m[11] +
    m[12] * m[7]  * m[10];

  inv[8] = m[4]  * m[9] * m[15] -
    m[4]  * m[11] * m[13] -
    m[8]  * m[5] * m[15] +
    m[8]  * m[7] * m[13] +
    m[12] * m[5] * m[11] -
    m[12] * m[7] * m[9];

  inv[12] = -m[4]  * m[9] * m[14] +
    m[4]  * m[10] * m[13] +
    m[8]  * m[5] * m[14] -
    m[8]  * m[6] * m[13] -
    m[12] * m[5] * m[10] +
    m[12] * m[6] * m[9];

  inv[1] = -m[1]  * m[10] * m[15] +
    m[1]  * m[11] * m[14] +
    m[9]  * m[2] * m[15] -
    m[9]  * m[3] * m[14] -
    m[13] * m[2] * m[11] +
    m[13] * m[3] * m[10];

  inv[5] = m[0]  * m[10] * m[15] -
    m[0]  * m[11] * m[14] -
    m[8]  * m[2] * m[15] +
    m[8]  * m[3] * m[14] +
    m[12] * m[2] * m[11] -
    m[12] * m[3] * m[10];

  inv[9] = -m[0]  * m[9] * m[15] +
    m[0]  * m[11] * m[13] +
    m[8]  * m[1] * m[15] -
    m[8]  * m[3] * m[13] -
    m[12] * m[1] * m[11] +
    m[12] * m[3] * m[9];

  inv[13] = m[0]  * m[9] * m[14] -
    m[0]  * m[10] * m[13] -
    m[8]  * m[1] * m[14] +
    m[8]  * m[2] * m[13] +
    m[12] * m[1] * m[10] -
    m[12] * m[2] * m[9];

  inv[2] = m[1]  * m[6] * m[15] -
    m[1]  * m[7] * m[14] -
    m[5]  * m[2] * m[15] +
    m[5]  * m[3] * m[14] +
    m[13] * m[2] * m[7] -
    m[13] * m[3] * m[6];

  inv[6] = -m[0]  * m[6] * m[15] +
    m[0]  * m[7] * m[14] +
    m[4]  * m[2] * m[15] -
    m[4]  * m[3] * m[14] -
    m[12] * m[2] * m[7] +
    m[12] * m[3] * m[6];

  inv[10] = m[0]  * m[5] * m[15] -
    m[0]  * m[7] * m[13] -
    m[4]  * m[1] * m[15] +
    m[4]  * m[3] * m[13] +
    m[12] * m[1] * m[7] -
    m[12] * m[3] * m[5];

  inv[14] = -m[0]  * m[5] * m[14] +
    m[0]  * m[6] * m[13] +
    m[4]  * m[1] * m[14] -
    m[4]  * m[2] * m[13] -
    m[12] * m[1] * m[6] +
    m[12] * m[2] * m[5];

  inv[3] = -m[1] * m[6] * m[11] +
    m[1] * m[7] * m[10] +
    m[5] * m[2] * m[11] -
    m[5] * m[3] * m[10] -
    m[9] * m[2] * m[7] +
    m[9] * m[3] * m[6];

  inv[7] = m[0] * m[6] * m[11] -
    m[0] * m[7] * m[10] -
    m[4] * m[2] * m[11] +
    m[4] * m[3] * m[10] +
    m[8] * m[2] * m[7] -
    m[8] * m[3] * m[6];

  inv[11] = -m[0] * m[5] * m[11] +
    m[0] * m[7] * m[9] +
    m[4] * m[1] * m[11] -
    m[4] * m[3] * m[9] -
    m[8] * m[1] * m[7] +
    m[8] * m[3] * m[5];

  inv[15] = m[0] * m[5] * m[10] -
    m[0] * m[6] * m[9] -
    m[4] * m[1] * m[10] +
    m[4] * m[2] * m[9] +
    m[8] * m[1] * m[6] -
    m[8] * m[2] * m[5];

  _T det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

  if (det == 0)
    return det;

  det = inverse(det);

  for (int i = 0; i < 16; i++)
    invOut[i] = inv[i] * det;

  return det;
}

template <typename _T>
Matrix4<_T> matrix4_inverse(const Matrix4<_T> &input) {
  Matrix4<_T> output;
  gluInvertMatrix((const _T *)input.v, (_T *)output.v);
  return output;
}

template <typename _T>
std::pair<Matrix4<_T>, _T> matrix4_inverse_det(const Matrix4<_T> &input) {
  Matrix4<_T> output;
  _T det = gluInvertMatrix((const _T *)input.v, (_T *)output.v);
  return std::make_pair(output, det);
}


template <typename _T>
Matrix3<_T> matrix3_inverse(const Matrix3<_T> &m) {
  Matrix3<_T> output;
  _T det = m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) -
           m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
           m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

  _T invdet = 1 / det;
  output[0][0] = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * invdet;
  output[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invdet;
  output[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invdet;
  output[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invdet;
  output[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invdet;
  output[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * invdet;
  output[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invdet;
  output[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * invdet;
  output[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invdet;
  return output;
}
#endif
