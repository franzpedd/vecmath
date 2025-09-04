#ifndef VECMATH_BASIC_OP_INCLUDED
#define VECMATH_BASIC_OP_INCLUDED

#include "vecmath_macros.h"
#include "vecmath_types.h"

/// @brief performs a+b, storing into result
void float2_add(const float2_t* a, const float2_t* b, float2_t* result);
void float3_add(const float3_t* a, const float3_t* b, float3_t* result);
void float4_add(const float4_t* a, const float4_t* b, float4_t* result);
void fmat2_add(const fmat2_t* a, const fmat2_t* b, fmat2_t* result);
void fmat3_add(const fmat3_t* a, const fmat3_t* b, fmat3_t* result);
void fmat4_add(const fmat4_t* a, const fmat4_t* b, fmat4_t* result);
//
void double2_add(const double2_t* a, const double2_t* b, double2_t* result);
void double3_add(const double3_t* a, const double3_t* b, double3_t* result);
void double4_add(const double4_t* a, const double4_t* b, double4_t* result);
void dmat2_add(const dmat2_t* a, const dmat2_t* b, dmat2_t* result);
void dmat3_add(const dmat3_t* a, const dmat3_t* b, dmat3_t* result);
void dmat4_add(const dmat4_t* a, const dmat4_t* b, dmat4_t* result);

/// @brief performs a-b, storing into result
void float2_sub(const float2_t* a, const float2_t* b, float2_t* result);
void float3_sub(const float3_t* a, const float3_t* b, float3_t* result);
void float4_sub(const float4_t* a, const float4_t* b, float4_t* result);
void fmat2_sub(const fmat2_t* a, const fmat2_t* b, fmat2_t* result);
void fmat3_sub(const fmat3_t* a, const fmat3_t* b, fmat3_t* result);
void fmat4_sub(const fmat4_t* a, const fmat4_t* b, fmat4_t* result);
//
void double2_sub(const double2_t* a, const double2_t* b, double2_t* result);
void double3_sub(const double3_t* a, const double3_t* b, double3_t* result);
void double4_sub(const double4_t* a, const double4_t* b, double4_t* result);
void dmat2_sub(const dmat2_t* a, const dmat2_t* b, dmat2_t* result);
void dmat3_sub(const dmat3_t* a, const dmat3_t* b, dmat3_t* result);
void dmat4_sub(const dmat4_t* a, const dmat4_t* b, dmat4_t* result);

/// @brief performs a * b, storing into result
void float2_mul(const float2_t* a, const float2_t* b, float2_t* result);
void float3_mul(const float3_t* a, const float3_t* b, float3_t* result);
void float4_mul(const float4_t* a, const float4_t* b, float4_t* result);
void fmat2_mul(const fmat2_t* a, const fmat2_t* b, fmat2_t* result);
void fmat3_mul(const fmat3_t* a, const fmat3_t* b, fmat3_t* result);
void fmat4_mul(const fmat4_t* a, const fmat4_t* b, fmat4_t* result);
//
void double2_mul(const double2_t* a, const double2_t* b, double2_t* result);
void double3_mul(const double3_t* a, const double3_t* b, double3_t* result);
void double4_mul(const double4_t* a, const double4_t* b, double4_t* result);
void dmat2_mul(const dmat2_t* a, const dmat2_t* b, dmat2_t* result);
void dmat3_mul(const dmat3_t* a, const dmat3_t* b, dmat3_t* result);
void dmat4_mul(const dmat4_t* a, const dmat4_t* b, dmat4_t* result);

/// @brief performs a / b, storing into result (let IEEE754 handle division by 0)
void float2_div(const float2_t* a, const float2_t* b, float2_t* result);
void float3_div(const float3_t* a, const float3_t* b, float3_t* result);
void float4_div(const float4_t* a, const float4_t* b, float4_t* result);
void fmat2_div(const fmat2_t* a, const fmat2_t* b, fmat2_t* result);
void fmat3_div(const fmat3_t* a, const fmat3_t* b, fmat3_t* result);
void fmat4_div(const fmat4_t* a, const fmat4_t* b, fmat4_t* result);
//
void double2_div(const double2_t* a, const double2_t* b, double2_t* result);
void double3_div(const double3_t* a, const double3_t* b, double3_t* result);
void double4_div(const double4_t* a, const double4_t* b, double4_t* result);
void dmat2_div(const dmat2_t* a, const dmat2_t* b, dmat2_t* result);
void dmat3_div(const dmat3_t* a, const dmat3_t* b, dmat3_t* result);
void dmat4_div(const dmat4_t* a, const dmat4_t* b, dmat4_t* result);

/// @brief checks if are exactly equals in value, this compare memory and is fast
vecbool float2_equals(const float2_t* a, const float2_t* b);
vecbool float3_equals(const float3_t* a, const float3_t* b);
vecbool float4_equals(const float4_t* a, const float4_t* b);
vecbool fmat2_equals(const fmat2_t* a, const fmat2_t* b);
vecbool fmat3_equals(const fmat3_t* a, const fmat3_t* b);
vecbool fmat4_equals(const fmat4_t* a, const fmat4_t* b);
//
vecbool double2_equals(const double2_t* a, const double2_t* b);
vecbool double3_equals(const double3_t* a, const double3_t* b);
vecbool double4_equals(const double4_t* a, const double4_t* b);
vecbool dmat2_equals(const dmat2_t* a, const dmat2_t* b);
vecbool dmat3_equals(const dmat3_t* a, const dmat3_t* b);
vecbool dmat4_equals(const dmat4_t* a, const dmat4_t* b);

/// @brief checks if almos-to-exactly equals in value, this compare value and is slow
vecbool float2_aprox_equals(const float2_t* a, const float2_t* b);
vecbool float3_aprox_equals(const float3_t* a, const float3_t* b);
vecbool float4_aprox_equals(const float4_t* a, const float4_t* b);
vecbool fmat2_aprox_equals(const fmat2_t* a, const fmat2_t* b);
vecbool fmat3_aprox_equals(const fmat3_t* a, const fmat3_t* b);
vecbool fmat4_aprox_equals(const fmat4_t* a, const fmat4_t* b);
//
vecbool double2_aprox_equals(const double2_t* a, const double2_t* b);
vecbool double3_aprox_equals(const double3_t* a, const double3_t* b);
vecbool double4_aprox_equals(const double4_t* a, const double4_t* b);
vecbool dmat2_aprox_equals(const dmat2_t* a, const dmat2_t* b);
vecbool dmat3_aprox_equals(const dmat3_t* a, const dmat3_t* b);
vecbool dmat4_aprox_equals(const dmat4_t* a, const dmat4_t* b);

#endif // VECMATH_BASIC_OP_INCLUDED