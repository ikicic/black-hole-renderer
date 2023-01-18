#ifndef FIELD_H
#define FIELD_H

class FlatDipole {
  double m;
  const double theta;
  const double phi;
  CartesianVector<double, 3> dipole;
 public:
  FlatDipole(
      double _surface_B = 2e11 * UNIT_T,
      double _theta = 0 * (M_PI / 180),
      double _phi = 0 * M_PI / 180)
        : m(1.6e30 * UNIT_J / UNIT_T), theta(_theta), phi(_phi) {
    _recalc_dipole();
    std::pair<double, double> surface_B = _surface_magnetic_field();
    m *= _surface_B / surface_B.second;
    _recalc_dipole();
    _debug_info();
  }

  inline void _recalc_dipole(void) {
    dipole = {{
      m * std::sin(theta) * std::cos(phi),
      m * std::sin(theta) * std::sin(phi),
      m * std::cos(theta)
    }};
  }

  template <typename T>
  inline CartesianVector4<T> get_potential_l(
      const CartesianVector4<T> &position_u) const {
    const CartesianVector<T, 3> position3{
      position_u[1],
      position_u[2],
      position_u[3],
    };
    const T r = position_u.get_r();
    return T(PHY_mu0 / (4 * M_PI)) * extend_vector(
        T(),
        -inverse(cube(r)) * cross(dipole, position3)
    );
  }

  template <typename Vector3>
  inline Vector3 _get_magnetic_field(const Vector3 &position) const {
    auto inv_r = inv_sqrt(position.sqr_length());
    return (PHY_mu0 / (4 * M_PI)) * cube(inv_r) * (
        3 * position.dot(dipole) * sqr(inv_r) * position
        - Vector3{{dipole[0], dipole[1], dipole[2]}}
    );
  }

  template <typename T>
  inline Matrix4<T> get_F_ll(const CartesianVector4<T> &position_u) const {
    const CartesianVector<T, 3> position3{
      position_u[1],
      position_u[2],
      position_u[3],
    };
    const CartesianVector<T, 3> B = _get_magnetic_field(position3);
    return Matrix4<T>{{
        {0, 0, 0, 0},
        {0, 0, -B[2], B[1]},
        {0, B[2], 0, -B[0]},
        {0, -B[1], B[0], 0}
    }};
  }

  inline double __magnetic_field(double r, double theta, double phi) const {
    CartesianVector<double, 3> position{
      r * std::sin(theta) * std::cos(phi),
      r * std::sin(theta) * std::sin(phi),
      r * std::cos(theta)
    };
    return _get_magnetic_field(position).length();
  }

  inline std::pair<double, double> _surface_magnetic_field(void) const {
    constexpr int N = 100;
    double sum = 0;
    double max = 0;
    for (int i = 0; i <= N; ++i) {
      double theta = random_double(0, M_PI);
      double phi = random_double(0, 2 * M_PI);
      double B = __magnetic_field(NEUTRON_STAR_r, theta, phi);
      max = std::max(B, max);
      sum += B;
    }

    fprintf(stderr, "_surface_magnetic_field  sum=%lg N=%d max=%lg  r=%lg\n", sum, N, max, NEUTRON_STAR_r);
    return std::make_pair(sum / N, max);
  }

  inline void _debug_info(void) const {
    constexpr int N = 21;
    for (int i = 0; i < N; ++i) {
      double theta = (1e-5 + (1 - 2e-5) * i) * M_PI / (N - 1);
      fprintf(stderr, "   theta=%.2lf deg B=%lg T\n",
          theta * 180 / M_PI,
          __magnetic_field(NEUTRON_STAR_r, theta, 0) / UNIT_T);
    }
    for (int i = 0; i < N; ++i) {
      double r = NEUTRON_STAR_r * (1 + 3. * i / (N / 1));
      fprintf(stderr, "   r=%.2lf deg B=%lg T\n",
          r / NEUTRON_STAR_r, __magnetic_field(r, M_PI / 2, 0) / UNIT_T);
    }

    std::pair<double, double> mag = _surface_magnetic_field();
    fprintf(stderr, "Neutron star r=%lg=%lg km\n",
        NEUTRON_STAR_r, NEUTRON_STAR_r / UNIT_km);
    fprintf(stderr, "Surface magnetic field  avg=%lg T=%lg\n",
        mag.first / UNIT_T, mag.first);
    fprintf(stderr, "Surface magnetic field  max=%lg T=%lg\n",
        mag.second / UNIT_T, mag.second);
    fprintf(stderr, "Magnetic dipole m=%lg J/T=%lg\n",
        m / UNIT_J * UNIT_T, m);
    std::cerr << "Dipole vector=" << dipole << '\n';
  }
};

#endif
