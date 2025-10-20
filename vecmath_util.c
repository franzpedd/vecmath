#include "vecmath_util.h"
#include <math.h>

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// angle utilities
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API float to_fradians(float degrees)
{
    return (float)(degrees * (VECMATH_EPSILON_PI / 180.0f));
}

VECMATH_API float to_fdegrees(float radians)
{
    return radians * (180.0f / (float)VECMATH_EPSILON_PI);
}

VECMATH_API float f_cos(float degree)
{
    return cosf(degree);
}

VECMATH_API float f_sin(float degree)
{
    return sinf(degree);
}

VECMATH_API float f_tan(float degree)
{
    return tanf(degree);
}

VECMATH_API double to_dradians(double degrees)
{
    return (double)(degrees * (VECMATH_EPSILON_PI / 180.0));
}

VECMATH_API double to_ddegrees(double radians)
{
    return radians * (180.0f / (float)VECMATH_EPSILON_PI);
}

VECMATH_API double d_cos(double degree)
{
    return cos(degree);
}

VECMATH_API double d_sin(double degree)
{
    return sin(degree);
}

VECMATH_API double d_tan(double degree)
{
    return tan(degree);
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// clamp utilities
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API float f_max(const float x, const float y)
{
    return (x < y) ? y : x;
}

VECMATH_API float f_min(const float x, const float y)
{
    return (y < x) ? y : x;
}

VECMATH_API float f_clamp(const float x, const float upper, const float lower)
{
    return f_min(upper, f_max(x, lower));
}

VECMATH_API double d_max(const double x, const double y)
{
    return (x < y) ? y : x;
}

VECMATH_API double d_min(const double x, const double y)
{
    return (y < x) ? y : x;
}

VECMATH_API double d_clamp(const double x, const double upper, const double lower)
{
    return f_min(upper, f_max(x, lower));
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// power
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API float f_power(float b, int e)
{
    // base case: anything^0 = 1
    if (e == 0) return 1.0f;
    
    // handle negative exponents
    if (e < 0) {
        if (e == VECMATH_EPSILON_INT_MIN) { // avoid overflow when negating INT_MIN
            return 1.0f / (b * f_power(b, -(e + 1)));
        }
        return 1.0f / f_power(b, -e);
    }
    
    // recursive case: use exponentiation by squaring
    float temp = f_power(b, e / 2);
    return (e % 2 == 0) ? temp * temp : b * temp * temp;
}

VECMATH_API double d_power(double b, int e)
{
    // base case: anything^0 = 1
    if (e == 0) return 1;
    
    // handle negative exponents
    if (e < 0) {
        if (e == VECMATH_EPSILON_INT_MIN) { // avoid overflow
            return 1 / (b * d_power(b, -(e + 1)));
        }
        return 1 / d_power(b, -e);
    }
    
    // recursive case: use exponentiation by squaring
    double temp = d_power(b, e / 2);
    return (e % 2 == 0) ? temp * temp : b * temp * temp;
}

VECMATH_API float f_log2(const float x)
{
    return log2f(x);
}

VECMATH_API double d_log2(const double x)
{
    return log2(x);
}

VECMATH_API int i_floor(const double x)
{
    if (x >= 0) return (int)x;
    
    int truncated = (int)x;
    return (x == truncated) ? truncated : truncated - 1;
}