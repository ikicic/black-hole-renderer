#ifndef BHR_CONFIG_H
#define BHR_CONFIG_H

namespace bhr {

// Configurable parameters.
typedef double real_t;
typedef double colreal_t;  // color real

#define DISK_DUMMY      1
#define DISK_KERTAP     2
#define DISK_SHAKURA    3

#define CHECK_KERR          0
#define CHECK_LAMBDA_PRECISION  0
#define GENERATE_LAMBDAS    0

#define PREPROCESS_LAMBDAS  1
#define RENDER_DISK         DISK_KERTAP
#define DISK_POLARIZATION   (RENDER_DISK == DISK_KERTAP)
#define SKY_ENABLED         0
#define MAGNETIC_FIELD      0
#define MAGNETIC_FIELD_FULL 0

#define FAKE_SKY      0
#define DISK_RELIEF_TEXTURE 0

#define PREDEFINED_PARAMS 1


// Checks and automatic macros.

#if MAGNETIC_FIELD_FULL && !MAGNETIC_FIELD
#error MAGNETIC_FIELD_FULL requires MAGNETIC_FIELD
#endif

#ifndef __clang__
# define CMATH_CONSTEXPR  constexpr
#else
# define CMATH_CONSTEXPR
#endif

/* Params BEGIN. */
#if PREDEFINED_PARAMS
# define PARAMS_CONSTEXPR        constexpr
# define PARAMS_CMATH_CONSTEXPR  CMATH_CONSTEXPR
# define PARAMS_FUNC_CONST
# define PARAMS_FUNC_STATIC                  static
# define PARAMS_FUNC_STATIC_CONSTEXPR        static constexpr
# define PARAMS_FUNC_STATIC_CMATH_CONSTEXPR  static CMATH_CONSTEXPR
#else
/* blah... */
# define PARAMS_CONSTEXPR                    const
# define PARAMS_CMATH_CONSTEXPR              const
# define PARAMS_FUNC_CONST                   const
# define PARAMS_FUNC_STATIC                  inline
# define PARAMS_FUNC_STATIC_CONSTEXPR        inline
# define PARAMS_FUNC_STATIC_CMATH_CONSTEXPR  inline
#endif
/* Params END. */

}  // namespace bhr

#endif
