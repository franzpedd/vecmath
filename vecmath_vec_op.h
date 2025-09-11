#ifndef VECMATH_VEC_OP_INCLUDED
#define VECMATH_VEC_OP_INCLUDED

#include "vecmath_defines.h"
#include "vecmath_types.h"

/// @ performs v * scalar whe n v is constant
VECMATH_API float2 float2_scalar(const float2* v, const float value);
VECMATH_API float3 float3_scalar(const float3* v, const float value);
VECMATH_API float4 float4_scalar(const float4* v, const float value);
VECMATH_API double2 double2_scalar(const double2* v, const float value);
VECMATH_API double3 double3_scalar(const double3* v, const float value);
VECMATH_API double4 double4_scalar(const double4* v, const float value);

/// @brief returns the length of the vector, (slower since performs sqrt)
VECMATH_API float float2_length(const float2* v);
VECMATH_API float float3_length(const float3* v);
VECMATH_API float float4_length(const float4* v);
VECMATH_API double double2_length(const double2* v);
VECMATH_API double double3_length(const double3* v);
VECMATH_API double double4_length(const double4* v);

/// @brief returns the length squared of the vector, (faster since doesn't call sqrt)
VECMATH_API float float2_length_sqrt(const float2* v);
VECMATH_API float float3_length_sqrt(const float3* v);
VECMATH_API float float4_length_sqrt(const float4* v);
VECMATH_API double double2_length_sqrt(const double2* v);
VECMATH_API double double3_length_sqrt(const double3* v);
VECMATH_API double double4_length_sqrt(const double4* v);

/// @brief normalizes the vectors, storing into result
VECMATH_API float2 float2_normalize(const float2* v);
VECMATH_API float3 float3_normalize(const float3* v);
VECMATH_API float4 float4_normalize(const float4* v);
VECMATH_API double2 double2_normalize(const double2* v);
VECMATH_API double3 double3_normalize(const double3* v);
VECMATH_API double4 double4_normalize(const double4* v);

/// @brief calculates the dot product
VECMATH_API float float2_dot(const float2* a, const float2* b);
VECMATH_API float float3_dot(const float3* a, const float3* b);
VECMATH_API float float4_dot(const float4* a, const float4* b);
VECMATH_API double double2_dot(const double2* a, const double2* b);
VECMATH_API double double3_dot(const double3* a, const double3* b);
VECMATH_API double double4_dot(const double4* a, const double4* b);

/// @brief calculates the cross product
VECMATH_API float float2_cross(const float2* a, const float2* b);
VECMATH_API float3 float3_cross(const float3* a, const float3* b);
VECMATH_API double double2_cross(const double2* a, const double2* b);
VECMATH_API double3 double3_cross(const double3* a, const double3* b);

/// @brief scales the vector v with scalar, storing into result
VECMATH_API float2 float2_scale(const float2* v, float scalar);
VECMATH_API float3 float3_scale(const float3* v, float scalar);
VECMATH_API float4 float4_scale(const float4* v, float scalar);
VECMATH_API double2 double2_scale(const double2* v, double scalar);
VECMATH_API double3 double3_scale(const double3* v, double scalar);
VECMATH_API double4 double4_scale(const double4* v, double scalar);

/// @brief performs linear interpolation between a and b, storing into result
VECMATH_API float2 float2_lerp(const float2* a, const float2* b, float t);
VECMATH_API float3 float3_lerp(const float3* a, const float3* b, float t);
VECMATH_API float4 float4_lerp(const float4* a, const float4* b, float t);
VECMATH_API double2 double2_lerp(const double2* a, const double2* b, double t);
VECMATH_API double3 double3_lerp(const double3* a, const double3* b, double t);
VECMATH_API double4 double4_lerp(const double4* a, const double4* b, double t);

/// @brief calculates the distance of two vectors, (slower since performs sqrt)
VECMATH_API float float2_distance(const float2* a, const float2* b);
VECMATH_API float float3_distance(const float3* a, const float3* b);
VECMATH_API float float4_distance(const float4* a, const float4* b);
VECMATH_API double double2_distance(const double2* a, const double2* b);
VECMATH_API double double3_distance(const double3* a, const double3* b);
VECMATH_API double double4_distance(const double4* a, const double4* b);

/// @brief calculates the distance squared of two vectors, (faster since doesn't call sqrt)
VECMATH_API float float2_distance_sqrt(const float2* a, const float2* b);
VECMATH_API float float3_distance_sqrt(const float3* a, const float3* b);
VECMATH_API float float4_distance_sqrt(const float4* a, const float4* b);
VECMATH_API double double2_distance_sqrt(const double2* a, const double2* b);
VECMATH_API double double3_distance_sqrt(const double3* a, const double3* b);
VECMATH_API double double4_distance_sqrt(const double4* a, const double4* b);

/// @brief calculates the reflection of a vector given it's normal, storing into result
VECMATH_API float2 float2_reflect(const float2* v, const float2* normal);
VECMATH_API float3 float3_reflect(const float3* v, const float3* normal);
VECMATH_API float4 float4_reflect(const float4* v, const float4* normal);
VECMATH_API double2 double2_reflect(const double2* v, const double2* normal);
VECMATH_API double3 double3_reflect(const double3* v, const double3* normal);
VECMATH_API double4 double4_reflect(const double4* v, const double4* normal);

/// @brief calculates the projection of two vectors, storing into result
VECMATH_API float2 float2_project(const float2* a, const float2* b);
VECMATH_API float3 float3_project(const float3* a, const float3* b);
VECMATH_API float4 float4_project(const float4* a, const float4* b);
VECMATH_API double2 double2_project(const double2* a, const double2* b);
VECMATH_API double3 double3_project(const double3* a, const double3* b);
VECMATH_API double4 double4_project(const double4* a, const double4* b);

#endif // VECMATH_VEC_OP_INCLUDED
