#ifndef VECMATH_BASIC_OP_INCLUDED
#define VECMATH_BASIC_OP_INCLUDED

#include "vecmath_macros.h"
#include "vecmath_types.h"

/// @brief performs a+b, storing into result
void float2_add(const float2* a, const float2* b, float2* result);
void float3_add(const float3* a, const float3* b, float3* result);
void float4_add(const float4* a, const float4* b, float4* result);
void fmat2_add(const fmat2* a, const fmat2* b, fmat2* result);
void fmat3_add(const fmat3* a, const fmat3* b, fmat3* result);
void fmat4_add(const fmat4* a, const fmat4* b, fmat4* result);
//
void double2_add(const double2* a, const double2* b, double2* result);
void double3_add(const double3* a, const double3* b, double3* result);
void double4_add(const double4* a, const double4* b, double4* result);
void dmat2_add(const dmat2* a, const dmat2* b, dmat2* result);
void dmat3_add(const dmat3* a, const dmat3* b, dmat3* result);
void dmat4_add(const dmat4* a, const dmat4* b, dmat4* result);

/// @brief performs a-b, storing into result
void float2_sub(const float2* a, const float2* b, float2* result);
void float3_sub(const float3* a, const float3* b, float3* result);
void float4_sub(const float4* a, const float4* b, float4* result);
void fmat2_sub(const fmat2* a, const fmat2* b, fmat2* result);
void fmat3_sub(const fmat3* a, const fmat3* b, fmat3* result);
void fmat4_sub(const fmat4* a, const fmat4* b, fmat4* result);
//
void double2_sub(const double2* a, const double2* b, double2* result);
void double3_sub(const double3* a, const double3* b, double3* result);
void double4_sub(const double4* a, const double4* b, double4* result);
void dmat2_sub(const dmat2* a, const dmat2* b, dmat2* result);
void dmat3_sub(const dmat3* a, const dmat3* b, dmat3* result);
void dmat4_sub(const dmat4* a, const dmat4* b, dmat4* result);

/// @brief performs a * b, storing into result
void float2_mul(const float2* a, const float2* b, float2* result);
void float3_mul(const float3* a, const float3* b, float3* result);
void float4_mul(const float4* a, const float4* b, float4* result);
void fmat2_mul(const fmat2* a, const fmat2* b, fmat2* result);
void fmat3_mul(const fmat3* a, const fmat3* b, fmat3* result);
void fmat4_mul(const fmat4* a, const fmat4* b, fmat4* result);
//
void double2_mul(const double2* a, const double2* b, double2* result);
void double3_mul(const double3* a, const double3* b, double3* result);
void double4_mul(const double4* a, const double4* b, double4* result);
void dmat2_mul(const dmat2* a, const dmat2* b, dmat2* result);
void dmat3_mul(const dmat3* a, const dmat3* b, dmat3* result);
void dmat4_mul(const dmat4* a, const dmat4* b, dmat4* result);

/// @brief performs a / b, storing into result (let IEEE754 handle division by 0)
void float2_div(const float2* a, const float2* b, float2* result);
void float3_div(const float3* a, const float3* b, float3* result);
void float4_div(const float4* a, const float4* b, float4* result);
void fmat2_div(const fmat2* a, const fmat2* b, fmat2* result);
void fmat3_div(const fmat3* a, const fmat3* b, fmat3* result);
void fmat4_div(const fmat4* a, const fmat4* b, fmat4* result);
//
void double2_div(const double2* a, const double2* b, double2* result);
void double3_div(const double3* a, const double3* b, double3* result);
void double4_div(const double4* a, const double4* b, double4* result);
void dmat2_div(const dmat2* a, const dmat2* b, dmat2* result);
void dmat3_div(const dmat3* a, const dmat3* b, dmat3* result);
void dmat4_div(const dmat4* a, const dmat4* b, dmat4* result);

/// @brief checks if are exactly equals in value, this compare memory and is fast
vecbool float2_equals(const float2* a, const float2* b);
vecbool float3_equals(const float3* a, const float3* b);
vecbool float4_equals(const float4* a, const float4* b);
vecbool fmat2_equals(const fmat2* a, const fmat2* b);
vecbool fmat3_equals(const fmat3* a, const fmat3* b);
vecbool fmat4_equals(const fmat4* a, const fmat4* b);
//
vecbool double2_equals(const double2* a, const double2* b);
vecbool double3_equals(const double3* a, const double3* b);
vecbool double4_equals(const double4* a, const double4* b);
vecbool dmat2_equals(const dmat2* a, const dmat2* b);
vecbool dmat3_equals(const dmat3* a, const dmat3* b);
vecbool dmat4_equals(const dmat4* a, const dmat4* b);

/// @brief checks if almos-to-exactly equals in value, this compare value and is slow
vecbool float2_aprox_equals(const float2* a, const float2* b);
vecbool float3_aprox_equals(const float3* a, const float3* b);
vecbool float4_aprox_equals(const float4* a, const float4* b);
vecbool fmat2_aprox_equals(const fmat2* a, const fmat2* b);
vecbool fmat3_aprox_equals(const fmat3* a, const fmat3* b);
vecbool fmat4_aprox_equals(const fmat4* a, const fmat4* b);
//
vecbool double2_aprox_equals(const double2* a, const double2* b);
vecbool double3_aprox_equals(const double3* a, const double3* b);
vecbool double4_aprox_equals(const double4* a, const double4* b);
vecbool dmat2_aprox_equals(const dmat2* a, const dmat2* b);
vecbool dmat3_aprox_equals(const dmat3* a, const dmat3* b);
vecbool dmat4_aprox_equals(const dmat4* a, const dmat4* b);

#endif // VECMATH_BASIC_OP_INCLUDED