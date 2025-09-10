#ifndef VECMATH_UTIL_INCLUDED
#define VECMATH_UTIL_INCLUDED

/// @brief angle utilities
float to_fradians(float degrees);
float to_fdegrees(float radians);
float f_cos(float degree);
float f_sin(float degree);
float f_tan(float degree);
//
double to_dradians(double degrees);
double to_ddegrees(double radians);
double d_cos(double degree);
double d_sin(double degree);
double d_tan(double degree);

/// @brief clamp utilities
float f_max(const float x, const float y);
float f_min(const float x, const float y);
float f_clamp(const float x, const float upper, const float lower);
//
double d_max(const double x, const double y);
double d_min(const double x, const double y);
double d_clamp(const double x, const double upper, const double lower);

/// @brief calculates b^e smartly (OlogN)
float f_power(float b, int e);
double d_power(double b, int e);

#endif // VECMATH_UTIL_INCLUDED