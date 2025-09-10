#include "vecmath_util.h"

#include "vecmath_macros.h"
#include <math.h>

#define UTIL_INT_MIN (-2147483647 - 1)

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// angle utilities
/////////////////////////////////////////////////////////////////////////////////////

float to_fradians(float degrees)
{
    return (float)(degrees * (VECMATH_EPSILON_PI / 180.0f));
}

float to_fdegrees(float radians)
{
    return radians * (180.0f / (float)VECMATH_EPSILON_PI);
}

float f_cos(float degree)
{
    return cosf(degree);
}

float f_sin(float degree)
{
    return sinf(degree);
}

float f_tan(float degree)
{
    return tanf(degree);
}

double to_dradians(double degrees)
{
    return (double)(degrees * (VECMATH_EPSILON_PI / 180.0));
}

double to_ddegrees(double radians)
{
    return radians * (180.0f / (float)VECMATH_EPSILON_PI);
}

double d_cos(double degree)
{
    return cos(degree);
}

double d_sin(double degree)
{
    return sin(degree);
}

double d_tan(double degree)
{
    return tan(degree);
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// clamp utilities
/////////////////////////////////////////////////////////////////////////////////////

float f_max(const float x, const float y)
{
    return (x < y) ? y : x;
}

float f_min(const float x, const float y)
{
    return (y < x) ? y : x;
}

float f_clamp(const float x, const float upper, const float lower)
{
    return f_min(upper, f_max(x, lower));
}

double d_max(const double x, const double y)
{
    return (x < y) ? y : x;
}

double d_min(const double x, const double y)
{
    return (y < x) ? y : x;
}

double d_clamp(const double x, const double upper, const double lower)
{
    return f_min(upper, f_max(x, lower));
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// power
/////////////////////////////////////////////////////////////////////////////////////

float f_power(float b, int e)
{
    // base case: anything^0 = 1
    if (e == 0) return 1.0f;
    
    // handle negative exponents
    if (e < 0) {
        if (e == UTIL_INT_MIN) { // avoid overflow when negating INT_MIN
            return 1.0f / (b * f_power(b, -(e + 1)));
        }
        return 1.0f / f_power(b, -e);
    }
    
    // recursive case: use exponentiation by squaring
    float temp = f_power(b, e / 2);
    return (e % 2 == 0) ? temp * temp : b * temp * temp;
}

double d_power(double b, int e)
{
    // base case: anything^0 = 1
    if (e == 0) return 1;
    
    // handle negative exponents
    if (e < 0) {
        if (e == UTIL_INT_MIN) { // avoid overflow
            return 1 / (b * d_power(b, -(e + 1)));
        }
        return 1 / d_power(b, -e);
    }
    
    // recursive case: use exponentiation by squaring
    double temp = d_power(b, e / 2);
    return (e % 2 == 0) ? temp * temp : b * temp * temp;
}