#ifndef DISK_H
#define DISK_H

#if RENDER_DISK

#include <bhr/chandra1960.h>
#include <bhr/colors.h>
#include <bhr/matrix.h>
#include <bhr/mod.h>
#include <bhr/parameters.h>
#include <bhr/physical_constants.h>
#include <bhr/spectrum.h>
#include <bhr/texture.h>

#include <tuple>

namespace bhr {

#if DISK_RELIEF_TEXTURE
extern Image disk_relief_tex;

bool load_disk_relief_texture(void);
#endif

class ShakuraSunyaevDisk;
extern ShakuraSunyaevDisk *_shakura_sunyaev;
double shakura_sunyaev_height(double r);

class ShakuraSunyaevDisk {
  double M;
  double Mdot;
  double Rstar;
  static constexpr double alpha = 0.1;   // ?
 public:
  ShakuraSunyaevDisk(double M, double Mdot, double Rstar)
      : M(M), Mdot(Mdot), Rstar(Rstar) {
    _shakura_sunyaev = this;
  }

  struct TexCoord {
    double r;               // r coordinate.
    double redshift_inf;    // Redshift assuming the observer is at infinity.
    double intensity;
    double T;
    double T4;
    double cosphi;
    double sinphi;

    friend inline TexCoord operator*(real_t c, const TexCoord &tex) {
      return {c * tex.r, c * tex.redshift_inf,
              c * tex.intensity, c * tex.T, c * tex.T4,
              c * tex.cosphi, c * tex.sinphi};
    }
    friend inline TexCoord operator+(const TexCoord &A, const TexCoord &B) {
      return {A.r + B.r, A.redshift_inf + B.redshift_inf,
              A.intensity + B.intensity,
              A.T + B.T, A.T4 + B.T4,
              A.cosphi + B.cosphi, A.sinphi + B.sinphi};
    }
  };

  struct Relativistic {
    double T, H;

    Relativistic(double a, double r, double M, double Mdot) {
      if (r < INNER_RADIUS) {
        T = 0;
        H = 0;
        return;
      }
      using std::log;
      using std::pow;
      double as = a / M;
      double as2 = as * as;
      double rs = r / M;
      double irs = 1 / rs;
      double irs2 = irs * irs;
      double irs3 = irs2 * irs;
      double irs32 = std::sqrt(irs2);
      double A = 1 + as2 * irs2 + 2 * as2 * irs3;
      double B = 1 + as * irs32;
      double C = 1 - 3 * irs + 2 * as * irs32;
      double D = 1 - 2 * irs + as2 * irs2;
      double E = 1 + 4 * as2 * irs2 - 4 * as2 * irs3 + 3 * sqr(as2 * irs2);
      // double F = 1 - 2 * as * irs32 + as2 * irs2;
      // double G = 1 - 2 * irs + as * irs32;
      double x = std::sqrt(r / M);
      double x0 = std::sqrt(_KERR_ISCO / M);
      double _phi = std::acos(as);
      double x1 = 2 * std::cos((_phi - M_PI) / 3);
      double x2 = 2 * std::cos((_phi + M_PI) / 3);
      double x3 = -2 * std::cos(_phi / 3);
      double _I = (1 + as / cube(x))  / (1 - 3 / sqr(x) + 2 * as / cube(x)) / x;
      double _II = x - x0 - 1.5 * as * log(x / x0)
          - 3 * sqr(x1 - as) / (x1 * (x1 - x2) * (x1 - x3)) * log((x - x1) / (x0 - x1))
          - 3 * sqr(x2 - as) / (x2 * (x2 - x3) * (x2 - x1)) * log((x - x2) / (x0 - x2))
          - 3 * sqr(x3 - as) / (x3 * (x3 - x1) * (x3 - x2)) * log((x - x3) / (x0 - x3));
      double I = _I * _II;
      double Mstar = M / (3 * PHY_M_Sun);
      double Mdot0star = Mdot / (1e14 * UNIT_kg / UNIT_s);

      T = 8e7 * UNIT_K
          * pow(alpha, -1. / 5)
          * pow(Mstar, -1. / 2)
          * pow(Mdot0star, 3. / 10)
          * pow(rs, -3. / 4)
          * pow(A, -1. / 10)
          * pow(B, -1. / 5)
          * pow(D, -3. / 20)
          * pow(E, 1. / 20)
          * pow(I, 3. / 10);
      H = 9e0 * UNIT_m
          * pow(alpha, -1. / 10)
          * pow(Mstar, 3. / 4)
          * pow(Mdot0star, 3. / 20)
          * pow(rs, 9. / 8)
          * pow(A, 19. / 20)
          * pow(B, -11. / 10)
          * pow(C, 1. / 2)
          * pow(D, -23. / 40)
          * pow(E, -19. / 40)
          * pow(I, 3. / 20);
    }
  };

  template <typename Spacetime, typename FullGeodesicData>
  TexCoord get_tex_coord(const Spacetime &spacetime,
                         const FullGeodesicData &geo) const {
    typedef decltype(geo.basic.position) Coord;
    typedef typename Coord::value_type U;

    Matrix4<U> tetrad_to_coord = spacetime.keplerian_tetrad_to_coord(
        geo.basic.position);
    Coord u = tetrad_to_coord * Coord{1, 0, 0, 0};
    double redshift_inf = -1 / spacetime.dot(
        geo.basic.position, geo.basic.direction, u);

    double r = geo.basic.position.get_r();
    Relativistic data(BLACK_HOLE_a, r, BLACK_HOLE_M, Mdot);
    double T = data.T;

    double intensity = PHY_black_body_sigma * sqr(sqr(T * redshift_inf));

    return TexCoord{
        r,
        redshift_inf,
        intensity,
        T,
        sqr(sqr(T)),
        std::cos(geo.basic.position.phi),
        std::sin(geo.basic.position.phi)
    };
  }

  double calc_total_luminosity(void) const {
    const int N = 10000;
    const double dr = (OUTER_RADIUS - INNER_RADIUS) / N;
    double luminosity = 0;
    for (int i = 0; i < N; ++i) {
      double r = INNER_RADIUS + (i + 1e-9) * dr;
      Relativistic data(BLACK_HOLE_a, r, BLACK_HOLE_M, Mdot);
      double T = data.T;
      double intensity = PHY_black_body_sigma * sqr(sqr(T));
      luminosity += intensity * 2 * M_PI * r * dr;
    }
    return 2 * luminosity;  // Two sides.
  }

  inline double get_height(double r) const {
    Relativistic data(BLACK_HOLE_a, r, BLACK_HOLE_M, Mdot);
    return data.H;
  }

  inline RGBd get_color(const TexCoord &tex) const {
    return discrete_linear_color_scale(21.7, log10(tex.intensity / UNIT_W), 27.8);
    // return discrete_linear_color_scale(4.5, log10(tex.T / UNIT_K), 6);

    // double c = 1e-27 * tex.intensity / UNIT_W;
    // return RGBd{c, c, c};

//    double T = std::pow(tex.T4, 0.25);
//#if DISK_RELIEF_TEXTURE
//    double phi = std::atan2(tex.sinphi, tex.cosphi);
//    double phix = Mod((phi / M_PI + 1 - .5) / 2, 1.);
//    // double x = (OUTER_RADIUS - tex.r) / (OUTER_RADIUS - INNER_RADIUS);
//    double x = tex.r / OUTER_RADIUS;
//    RGBA tex_color = disk_relief_tex.get_pixel_rel(0.3 / (.4 + x), phix);
//    T *= 1 + 0.1 * (tex_color.to_RGBd().luminance() - .5);
//#endif
//    RGBd result = get_black_body_color(T, tex.redshift_inf);
//    return 1e-15 * result;
  }

  inline RGBd get_color(const TexCoord &t1,
                        const TexCoord &t2,
                        const TexCoord &t3,
                        const TexCoord &t4) const {
    return 0.25 * (get_color(t1)
                 + get_color(t2)
                 + get_color(t3)
                 + get_color(t4));
  }
};



template <typename Coord>
class KERTAPDiskTexture {
  typedef typename Coord::value_type T;
  static constexpr real_t Gamma = 2.0;
  static constexpr real_t nr = 3.0;
  const Coord camera_position;
 public:

  KERTAPDiskTexture(const Coord &_camera_position)
      : camera_position(_camera_position) {}

  struct TexCoord {
    T r;               // r coordinate.
    T weight;
    T delta;
    T chi;
    T chi2_cos;
    T chi2_sin;
    T mu;
    T redshift_inf;    // Redshift assuming the observer is at infinity.
    T intensity;

    friend inline TexCoord operator*(real_t c, const TexCoord &tex) {
      return {c * tex.r, c * tex.weight, c * tex.delta, c * tex.chi,
              c * tex.chi2_cos, c * tex.chi2_sin,
              c * tex.mu, c * tex.redshift_inf, c * tex.intensity};
    }
    friend inline TexCoord operator+(const TexCoord &A, const TexCoord &B) {
      return {A.r + B.r, A.weight + B.weight, A.delta + B.delta, A.chi + B.chi,
              A.chi2_cos + B.chi2_cos, A.chi2_sin + B.chi2_sin,
              A.mu + B.mu, A.redshift_inf + B.redshift_inf,
              A.intensity + B.intensity};
    }
  };

  /* For forward parallel transport (used for testing transport matrix M). */
  struct _ForwardsState {
    Coord position;
    Coord direction;
    Coord vec;

    template <typename U>
    friend inline void mult(_ForwardsState *self, const U &c) {
      mult(&self->position, c);
      mult(&self->direction, c);
      mult(&self->vec, c);
    }

    template <typename U>
    friend inline void mult_add(
        _ForwardsState *self, const U &c, const _ForwardsState &A) {
      mult_add(&self->position, c, A.position);
      mult_add(&self->direction, c, A.direction);
      mult_add(&self->vec, c, A.vec);
    }

    template <typename U>
    friend inline void set_and_mult_add(
        _ForwardsState *self,
        const _ForwardsState &A, const U &c, const _ForwardsState &B) {
      set_and_mult_add(&self->position, A.position, c, B.position);
      set_and_mult_add(&self->direction, A.direction, c, B.direction);
      set_and_mult_add(&self->vec, A.vec, c, B.vec);
    }

    friend inline double numerical_distance(_ForwardsState &A,
                                            _ForwardsState &B) {
      return numerical_distance(A.position, B.position)
           + numerical_distance(A.direction, B.direction)
           + numerical_distance(A.vec, B.vec);
    }
  };

  template <typename Spacetime>
  _ForwardsState _simulate_forwards(const Spacetime &spacetime,
                                    const std::vector<double> &dlambdas,
                                    _ForwardsState state) const {
    auto RHS = [&spacetime](const _ForwardsState &state) {
      _ForwardsState result;
      auto christoffel_ull = spacetime.get_christoffel_ull(state.position);
      for (int k = 0; k < 4; ++k) {
        T tmp1 = T();
        T tmp2 = T();
        for (int i = 0; i < 4; ++i)
          for (int j = 0; j < 4; ++j) {
            tmp1 += state.direction[i]
                  * christoffel_ull[k][i][j]
                  * state.direction[j];
            tmp2 += state.direction[i]
                  * christoffel_ull[k][i][j]
                  * state.vec[j];
          }
        result.direction[k] = -tmp1;
        result.vec[k] = -tmp2;
      }
      result.position = state.direction;
      return result;
    };

    _ForwardsState new_state;
    for (int s = (int)dlambdas.size() - 1; s >= 0; --s) {
      double dlambda = dlambdas[s];
      integration_step__RGF45(RHS, dlambda, state, &new_state);
      state = new_state;
    }
    return state;
  }

  template <typename Spacetime, typename FullGeodesicData>
  TexCoord get_tex_coord(const Spacetime &spacetime,
                         const FullGeodesicData &geo) const {
    Matrix4<T> tetrad_to_coord = spacetime.keplerian_tetrad_to_coord(
        geo.basic.position);
    Matrix4<T> coord_to_tetrad = matrix4_inverse(tetrad_to_coord);

    /* Polarization in the source local coordinate system. */
    Coord photon_hat = coord_to_tetrad * geo.basic.direction;
    T E_source__norm = std::sqrt(sqr(photon_hat.phi) + sqr(photon_hat.r));
    Coord E_source_hat{0, -photon_hat.phi / E_source__norm,
                        0, photon_hat.r / E_source__norm};
    Coord E_source = tetrad_to_coord * E_source_hat;

    // Coord start_position = camera_position;
    Coord start_position = geo.extra.start_position;

    /* Parallel transport of polarization from source to the observer. */
    Matrix4<T> transport_mat = matrix4_inverse(
        geo.basic.parallel_transport_lu);
    Coord E_observer = E_source * transport_mat;  // Yes, this order.
    Coord E_observer_hat =
        spacetime.ZAMO_coord_to_tetrad(start_position) * E_observer;

    // T chi = std::atan(-E_observer_hat.theta / E_observer_hat.phi);
    T chi = std::atan2(-E_observer_hat.theta, E_observer_hat.phi);
    T chi2_cos = std::cos(2 * chi);
    T chi2_sin = std::sin(2 * chi);

    /* Chandra table values for weight and delta. */
    T norm = std::sqrt(
        sqr(photon_hat.r) + sqr(photon_hat.theta) + sqr(photon_hat.phi));
    T mu = std::abs(photon_hat.theta) / norm;
    T chandra_weight, chandra_delta;
    std::tie(chandra_weight, chandra_delta) = chandra1960(mu);

    /* Redshift calculation. */
    // u_hat = {1, 0, 0, 0}
    Coord u = tetrad_to_coord * Coord{1, 0, 0, 0};
    T redshift_inf = -1 / spacetime.dot(
        geo.basic.position, geo.basic.direction, u);

    // Redshifts goes like g^-4, not g^-3! (?)
    // (g = nu_obs / nu_src)
    T intensity = chandra_weight
                   * std::pow(redshift_inf, this->Gamma - 1 + 4)
                   / std::pow(geo.basic.position.get_r(), this->nr);

    return TexCoord{
        geo.basic.position.get_r(),
        chandra_weight,
        chandra_delta,
        chi,
        chi2_cos,
        chi2_sin,
        mu,
        redshift_inf,
        intensity
    };
  }

  inline RGBd get_color(const TexCoord &tex) const {
    // return linear_color_scale(-M_PI / 2, tex.chi, M_PI / 2);
    return discrete_linear_color_scale(0, log10(tex.intensity) + 2.3, 6.0);
    // return linear_color_scale(0, tex.chi, 0.2);
    // return linear_color_scale(0.0, tex.redshift_inf, 1.5);
    // return linear_color_scale(0.0, tex.delta, .12);
    // return linear_color_scale(0.0, tex.mu, 1.0);
  }

  static constexpr double _arrow_zoom = 30;
  inline std::pair<T, T> get_arrow_vector(const TexCoord &tex) const {
    return std::make_pair(_arrow_zoom * tex.delta * std::cos(tex.chi),
                          _arrow_zoom * tex.delta * std::sin(tex.chi));
  }
  inline T _get_arrow_scale_01(void) const {
    return _arrow_zoom * 0.1;
  }

  inline RGBd get_color(const TexCoord &t1,
                        const TexCoord &t2,
                        const TexCoord &t3,
                        const TexCoord &t4) const {
    return 0.25 * (get_color(t1)
                 + get_color(t2)
                 + get_color(t3)
                 + get_color(t4));
  }
};



class DummyDiskTexture {
 public:
  DummyDiskTexture() {}

  struct TexCoord {
    real_t r;   // radius
    real_t cosphi;
    real_t sinphi;

    friend inline TexCoord operator*(real_t c, const TexCoord &tex) {
      return {c * tex.r, c * tex.cosphi, c * tex.sinphi};
    }
    friend inline TexCoord operator+(const TexCoord &A, const TexCoord &B) {
      return {A.r + B.r, A.cosphi + B.cosphi, A.sinphi + B.sinphi};
    }
  };

  template <typename Spacetime, typename FullGeodesicData>
  TexCoord get_tex_coord(const Spacetime &spacetime,
                         const FullGeodesicData &geo) const {
    const auto &position = geo.basic.position;
    const auto &direction = geo.basic.direction;
    const real_t r = position.get_r();

    // Approximate value of the frequency factor.
    const auto spherical = position.spherical_part(
        spacetime.coord_system_parameters(position));

    return TexCoord{r, std::cos(spherical.phi), std::sin(spherical.phi)};
  }

  inline RGBd get_color(const TexCoord &tex) const {
    double x = double((double(OUTER_RADIUS) - tex.r) / (OUTER_RADIUS - INNER_RADIUS));
    double phi = atan2(tex.sinphi, tex.cosphi);
    phi = phi - Mod(phi, 2 * M_PI / 24);
    return RGBd(
        .5 + .5 * cos(phi),
        .25 + (int)(x * 10) / 10. * .75,
        0
    );
  }

  inline RGBd get_color(const TexCoord &t1,
                        const TexCoord &t2,
                        const TexCoord &t3,
                        const TexCoord &t4) const {
    return 0.25 * (get_color(t1)
                 + get_color(t2)
                 + get_color(t3)
                 + get_color(t4));
  }
};

#endif

}  // namespace bhr

#endif
