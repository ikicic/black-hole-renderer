#ifndef QED_LAGRANGIAN_H
#define QED_LAGRANGIAN_H

#include <bhr/autodiff_2nd.h>
#include <bhr/numerics.h>
#include <tuple>

namespace QED {
  // REFERENCE: 1101.3433v5 - Kim, Lee - Light bending by nonlinear
  // electrodynamics under strong electric and magnetic field.
  // REFERENCE: 1501.06234 - Non-linear effects on radiation propagation around
  // a charged compact object
  // NOTE: 1501.06234 uses L = +a F + b (...).
  //
  // L = -a F + b (4 F^2 + 7 G^2)
  constexpr double _a = sqr(PHY_c) * PHY_eps0;
  constexpr double _b = 2 * sqr(PHY_alpha * PHY_eps0) * cube(PHY_hbar)
      / (45 * sqr(sqr(PHY_m_e)) * PHY_c);
  constexpr double lambda1 = 8 * _b / _a;
  constexpr double lambda2 = 14 * _b / _a;

  template <typename _T>
  _T lagrangian_lowest_order(const _T &F, const _T &G) {
    return _b * (4 * sqr(F) + 7 * sqr(G));
  };

#if GENERATE_LAMBDAS
  void generate_lambdas(void);
#endif
};

// WARNING: We can't calculate derivatives here!
// 1503.00532 results might not be correct, as the functions are used outside
// of their domains. (?)

template <typename _T>
_T magnetic_only_lagrangian__base_low(const _T &x) {
  static const _T coef[] = {
    0.02222222222222222222222222222222222222222,  // x^4
    -0.0126984126984126984126984126984126984127,  // x^6
    0.0253968253968253968253968253968253968254,   // x^8
    -0.1077441077441077441077441077441077441077,
    0.7854190254190254190254190254190254190254,
    -8.752136752136752136752136752136752136752,
    138.3308309990662931839402427637721755369,
    -2943.29339689092011073435222042033187544,
    81115.05579301789828105617579301789828106,
    -2.810784313595583160800552104899930986888e6,
    1.196126159272688663993011819098775620515e8,
    -6.132361388942222222222222222222222222222e9,
    3.72802864503590739103382781543701083931e11,
    -2.651652482863668848114295436717175156711e13
  };
  constexpr int N = sizeof(coef) / sizeof(coef[0]);
  _T xx = sqr(x);
  return sqr(xx) * evaluate_polynomial(coef, N, xx);
}

template <typename _T>
_T magnetic_only_lagrangian__base_middle(const _T &x) {
  static const _T coef[] = {
    -2.402822541120807862746802014035460672084e-11, // x^0
    1.581091978782811783512679891854091345588e-9,   // x^1
    -5.045130989201587775864983624930460474351e-8,  // x^2
    1.039016427150383697465477707768220928832e-6,
    0.02220672183332068779727928023653333213048,
    0.0001781189572732986671273146190196179506034,
    -0.01433378727010454721501994631413918466174,
    0.01226735497395546088889555300496929727068,
    -0.05071683889625349166031460117582256515388,
    0.3920964140076371776421692264288313590459,
    -1.776078427702768334700857351481566774042,
    5.768133944491510843986321190765196568575,
    -14.92526853311914887589146839925505691191,
    32.28664474951610653082866771467176810936,
    -59.74863232788376188498696352517538515546,
    95.58337096250938553765104129984057881506,
    -132.5948746697026941787044057702830887938,
    159.1344494423253396118551175694664745274,
    -163.9968892234770263241238278496753843607,
    143.0506287419047651994188139266341611458,
    -102.8427016220056981183840272163945044212,
    57.64064495842531832017234480952985551107,
    -21.44889422436460995797504963409553748671,
    0.9191261829268426489599456024726059179995,
    5.879427083102159680181614887195990050902,
    -5.247132184071285412228271986289316473751,
    2.753689963541177147012377065458977741492,
    -0.9878258188818254818533183630096398902582,
    0.2406987805075494911677265868961066425287,
    -0.03634343136003851511381405566268196424206,
    0.002586573473881441646672063552334009052919
  };
  constexpr int N = sizeof(coef) / sizeof(coef[0]);
  return evaluate_polynomial(coef, N, x);
}

template <typename _T>
_T magnetic_only_lagrangian__base_high(const _T &x) {
  _T y = inverse(x);
  static const _T coef[] = {
    0.8079657578292062244053600156878870685167,   // y^0
    0.1370778389040188697060345972205020991016,   // y^1
    -0.02504285214915821427916121169815520814094,
    0.006764520210694613696975023103382299392342,
    -0.002160266156548687346523678096785487850119,
    0.000756951683024143705144730602523006345165,
    -0.0002813474546266525744530685127929120422991,
    0.0001089493659068950021027219225812339914561,
    -0.00004348994760529870722299708199793455123635,
    0.00001777334117769563361749193819070168709173,
    -7.401856864081139504606876498331483313126e-6,
    3.130778265704214393963585480668255680841e-6,
    -1.341596617106448651521119286205729250346e-6,
    5.81322805138031706211204046401329597835e-7,
    -2.543209300324267121616495077796033382639e-7,
    1.121986930240640810654284052554452247646e-7,
    -4.986570456499276253505899268722848476899e-8,
    2.23082563008373104305123772346679941753e-8,
    -1.003869617073376860306313822854488876095e-8,
    4.541310600834954910260668945488173391006e-9,
    -2.06423110659088776288396479602418291217e-9,
    9.423661500102472189697729303414517739968e-10,
    -4.319177672568982776140272910131574032645e-10,
    1.986821610943851970154999916633383305473e-10,
    -9.169945623356594623531563337616783558094e-11,
    4.245345132734509055961915717583888916457e-11,
    -1.971053082655015031645359288094938137631e-11,
    9.175591902314937841316538774057011761417e-12,
    -4.281942879771106045791853632066741180732e-12,
    2.002844248350203823217437518203300043422e-12,
    -9.388332409769767062627552746552179029771e-13
  };
  constexpr int N = sizeof(coef) / sizeof(coef[0]);
  _T result = evaluate_polynomial(coef, N, y);
  result -= x * (0.1447298858494001741434273513530587116473
               + x * 0.7639688479484886137166012671517303816975);
  result -= std::log(y) * (0.5 + x * (1 + (_T(1.) / 3) * x));
  return result;
}

template <typename _T>
inline _T magnetic_only_lagrangian__base(const _T &x) {
  if (x < .13)
    return magnetic_only_lagrangian__base_low(x);
  else if (x < 0.97)
    return magnetic_only_lagrangian__base_middle(x);
  else
    return magnetic_only_lagrangian__base_high(x);
}

inline double _magnetic_only_lagrangian(double F) {
  /* NOTE: The formulas here already assumme G = 0, so no automatic
   * differentiation with respect to G is possible. */
  // REFERENCE: 0308223v1 (not actually used)
  // REFERENCE: PhysRevD.19.2929 (formula)
  // REFERENCE: 1101.3433v5 - Kim, Lee (factor)
  double B = std::sqrt(2 * F);
  double x = B / PHY_Bc;
  constexpr double factor = sqr(sqr(PHY_m_e * PHY_c)) * PHY_c
      / cube(PHY_hbar) / (8 * sqr(M_PI));
  return factor * magnetic_only_lagrangian__base(x);
}

// WARNING: No automatic differentiation!
inline double magnetic_only_lagrangian(double F, double /*G*/) {
  return _magnetic_only_lagrangian(F);
}

template <typename _T>
std::pair<_T, _T> qed_metric_correction_lambda__dimless(
    auto dimless_lagrangian_func, const _T &F, const _T &G) {
  // REFERENCE: http://arxiv.org/abs/1501.06234v2
  // NOTE: The article uses F = -F_ab F^ab / 4, while we use F = +F_ab F^ab/4.
  //                        G = -F_ab *F^ab / 4              G = +F_ab *F^ab/4.
  // So, we had to replace each occurence of F and d/dF with -F and -d/dF.

  using std::sqrt;
  typedef second_partial_derivatives<_T, 2> AD2;

  AD2 Fad{F, 1, 0, 0, 0, 0};
  AD2 Gad{G, 0, 1, 0, 0, 0};
  AD2 L = dimless_lagrangian_func(Fad, Gad);

  // Manually add -F to the lagrangian.
  // L.value() -= F;
  L.first(0) -= 1;

  _T det = L.second(0, 0) * L.second(1, 1) - sqr(L.second(0, 1));
  _T delta1 = -2 * F * det - L.first(0) * (L.second(1, 1) - L.second(0, 0));
  _T delta2 = -2 * G * det + 2 * L.first(0) * L.second(0, 1);
  _T delta = sqr(delta1) + sqr(delta2);
  _T delta_sqrt = sqrt(delta);

  _T num = L.first(0) * (L.second(0, 0) + L.second(1, 1)) - 2 * F * det;
  _T den = 2 * sqr(G) * det
         + 4 * L.first(0) * (L.second(1, 1) * F - L.second(0, 1) * G)
         - 2 * sqr(L.first(0));

  // den seems to be negative.
  _T lambda1 = (num + delta_sqrt) / den;
  _T lambda2 = (num - delta_sqrt) / den;

  return std::make_pair(lambda1, lambda2);
}

template <typename _T>
std::pair<_T, _T> qed_metric_correction_lambda(
    auto lagrangian_func, const _T &F, const _T &G) {
  std::pair<_T, _T> result = qed_metric_correction_lambda__dimless(
      [&lagrangian_func](auto dimlessF, auto dimlessG) {
        return PHY_inv_sqr_Bc * lagrangian_func(
            dimlessF * PHY_sqr_Bc, dimlessG * PHY_sqr_Bc);
      },
      F * PHY_inv_sqr_Bc,
      G * PHY_inv_sqr_Bc
  );
  result.first *= PHY_inv_sqr_Bc;
  result.second *= PHY_inv_sqr_Bc;
  return result;
}

namespace QEDCache {
  typedef first_partial_derivatives<real_t, 2> lambda_t;
  struct LambdaCacheType {
    lambda_t lambda[2];
    friend inline LambdaCacheType operator*(double A,
                                            const LambdaCacheType &B) {
      return {{A * B.lambda[0], A * B.lambda[1]}};
    }
    friend inline LambdaCacheType operator+(const LambdaCacheType &A,
                                            const LambdaCacheType &B) {
      return {{A.lambda[0] + B.lambda[0], A.lambda[1] + B.lambda[1]}};
    }
  };

  bool lambda_preprocess(int thread_count);
  LambdaCacheType get_cached_lambdas(double F, double G);
  void check_lambda_interpolation_precision(void);
}

#endif
