#ifndef VECMATH_DEFINES_INCLUDED
#define VECMATH_DEFINES_INCLUDED

/// @brief used macro definitions
#define VECMATH_EPSILON_FZERO 1e-6f
#define VECMATH_EPSILON_DZERO 1e-12

#define VECMATH_EPSILON_PI 3.14159265358979323846
#define VECMATH_EPSILON_FLT 1.192092896e-07F

#define VECMATH_EPSILON_INT_MIN (-2147483647 - 1)

/// @brief when using the header-only version this will be defined
#ifdef VECMATH_REQUESTING_HEADER_ONLY
    #define VECMATH_API static
#else
    #define VECMATH_API
#endif

#endif // VECMATH_DEFINES_INCLUDED