#ifndef VECMATH_VEC_OP_INCLUDED
#define VECMATH_VEC_OP_INCLUDED

#include "vecmath_types.h"

/// @brief returns the length of the vector, (slower since performs sqrt)
float float2_length(const float2* v);
float float3_length(const float3* v);
float float4_length(const float4* v);
double double2_length(const double2* v);
double double3_length(const double3* v);
double double4_length(const double4* v);

/// @brief returns the length squared of the vector, (faster since doesn't call sqrt)
float float2_length_sqrt(const float2* v);
float float3_length_sqrt(const float3* v);
float float4_length_sqrt(const float4* v);
double double2_length_sqrt(const double2* v);
double double3_length_sqrt(const double3* v);
double double4_length_sqrt(const double4* v);

/// @brief normalizes the vectors, storing into result
void float2_normalize(const float2* v, float2* result);
void float3_normalize(const float3* v, float3* result);
void float4_normalize(const float4* v, float4* result);
void double2_normalize(const double2* v, double2* result);
void double3_normalize(const double3* v, double3* result);
void double4_normalize(const double4* v, double4* result);

/// @brief calculates the dot product
float float2_dot(const float2* a, const float2* b);
float float3_dot(const float3* a, const float3* b);
float float4_dot(const float4* a, const float4* b);
double double2_dot(const double2* a, const double2* b);
double double3_dot(const double3* a, const double3* b);
double double4_dot(const double4* a, const double4* b);

/// @brief calculates the cross product
void float2_cross(const float2* a, const float2* b, float* result);
void float3_cross(const float3* a, const float3* b, float3* result);
void double2_cross(const double2* a, const double2* b, double* result);
void double3_cross(const double3* a, const double3* b, double3* result);

/// @brief scales the vector v with scalar, storing into result
void float2_scale(const float2* v, float scalar, float2* result);
void float3_scale(const float3* v, float scalar, float3* result);
void float4_scale(const float4* v, float scalar, float4* result);
void double2_scale(const double2* v, double scalar, double2* result);
void double3_scale(const double3* v, double scalar, double3* result);
void double4_scale(const double4* v, double scalar, double4* result);

/// @brief performs linear interpolation between a and b, storing into result
void float2_lerp(const float2* a, const float2* b, float t, float2* result);
void float3_lerp(const float3* a, const float3* b, float t, float3* result);
void float4_lerp(const float4* a, const float4* b, float t, float4* result);
void double2_lerp(const double2* a, const double2* b, double t, double2* result);
void double3_lerp(const double3* a, const double3* b, double t, double3* result);
void double4_lerp(const double4* a, const double4* b, double t, double4* result);

/// @brief calculates the distance of two vectors, (slower since performs sqrt)
float float2_distance(const float2* a, const float2* b);
float float3_distance(const float3* a, const float3* b);
float float4_distance(const float4* a, const float4* b);
double double2_distance(const double2* a, const double2* b);
double double3_distance(const double3* a, const double3* b);
double double4_distance(const double4* a, const double4* b);

/// @brief calculates the distance squared of two vectors, (faster since doesn't call sqrt)
float float2_distance_sqrt(const float2* a, const float2* b);
float float3_distance_sqrt(const float3* a, const float3* b);
float float4_distance_sqrt(const float4* a, const float4* b);
double double2_distance_sqrt(const double2* a, const double2* b);
double double3_distance_sqrt(const double3* a, const double3* b);
double double4_distance_sqrt(const double4* a, const double4* b);

/// @brief calculates the reflection of a vector given it's normal, storing into result
void float2_reflect(const float2* v, const float2* normal, float2* result);
void float3_reflect(const float3* v, const float3* normal, float3* result);
void float4_reflect(const float4* v, const float4* normal, float4* result);
void double2_reflect(const double2* v, const double2* normal, double2* result);
void double3_reflect(const double3* v, const double3* normal, double3* result);
void double4_reflect(const double4* v, const double4* normal, double4* result);

/// @brief calculates the projection of two vectors, storing into result
void float2_project(const float2* a, const float2* b, float2* result);
void float3_project(const float3* a, const float3* b, float3* result);
void float4_project(const float4* a, const float4* b, float4* result);
void double2_project(const double2* a, const double2* b, double2* result);
void double3_project(const double3* a, const double3* b, double3* result);
void double4_project(const double4* a, const double4* b, double4* result);

#endif // VECMATH_VEC_OP_INCLUDED