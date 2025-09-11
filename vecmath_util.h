#ifndef VECMATH_UTIL_INCLUDED
#define VECMATH_UTIL_INCLUDED

#include "vecmath_defines.h"

/// @brief angle utilities
VECMATH_API float to_fradians(float degrees);
VECMATH_API float to_fdegrees(float radians);
VECMATH_API float f_cos(float degree);
VECMATH_API float f_sin(float degree);
VECMATH_API float f_tan(float degree);
//
VECMATH_API double to_dradians(double degrees);
VECMATH_API double to_ddegrees(double radians);
VECMATH_API double d_cos(double degree);
VECMATH_API double d_sin(double degree);
VECMATH_API double d_tan(double degree);

/// @brief clamp utilities
VECMATH_API float f_max(const float x, const float y);
VECMATH_API float f_min(const float x, const float y);
VECMATH_API float f_clamp(const float x, const float upper, const float lower);
//
VECMATH_API double d_max(const double x, const double y);
VECMATH_API double d_min(const double x, const double y);
VECMATH_API double d_clamp(const double x, const double upper, const double lower);

/// @brief calculates b^e smartly (OlogN)
VECMATH_API float f_power(float b, int e);
VECMATH_API double d_power(double b, int e);

#endif // VECMATH_UTIL_INCLUDED