#ifndef VECMATH_VEC_OP_INCLUDED
#define VECMATH_VEC_OP_INCLUDED

#include "vecmath_types.h"

/// @brief returns the length of the vector, (slower since performs sqrt)
float float2_length(const float2_t* v);
float float3_length(const float3_t* v);
float float4_length(const float4_t* v);

/// @brief returns the length squared of the vector, (faster since doesn't call sqrt)
float float2_length_sqrt(const float2_t* v);
float float3_length_sqrt(const float3_t* v);
float float4_length_sqrt(const float4_t* v);

/// @brief normalizes the vectors, storing into result
void float2_normalize(const float2_t* v, float2_t* result);
void float3_normalize(const float3_t* v, float3_t* result);
void float4_normalize(const float4_t* v, float4_t* result);

/// @brief calculates the dot product
float float2_dot(const float2_t* a, const float2_t* b);
float float3_dot(const float3_t* a, const float3_t* b);
float float4_dot(const float4_t* a, const float4_t* b);

/// @brief calculates the cross product, 3D only, storing into result
void float3_cross(const float3_t* a, const float3_t* b, float3_t* result);

/// @brief scales the vector v with scalar, storing into result
void float2_scale(const float2_t* v, float scalar, float2_t* result);
void float3_scale(const float3_t* v, float scalar, float3_t* result);
void float4_scale(const float4_t* v, float scalar, float4_t* result);

/// @brief performs linear interpolation between a and b, storing into result
void float2_lerp(const float2_t* a, const float2_t* b, float t, float2_t* result);
void float3_lerp(const float3_t* a, const float3_t* b, float t, float3_t* result);
void float4_lerp(const float4_t* a, const float4_t* b, float t, float4_t* result);

/// @brief calculates the distance of two vectors
float float2_distance(const float2_t* a, const float2_t* b);
float float3_distance(const float3_t* a, const float3_t* b);
float float4_distance(const float4_t* a, const float4_t* b);

/// @brief calculates the distance squared of two vectors
float float2_distance_sqrt(const float2_t* a, const float2_t* b);
float float3_distance_sqrt(const float3_t* a, const float3_t* b);
float float4_distance_sqrt(const float4_t* a, const float4_t* b);

/// @brief calculates the reflection of a vector given it's normal, storing into result
void float2_reflect(const float2_t* v, const float2_t* normal, float2_t* result);
void float3_reflect(const float3_t* v, const float3_t* normal, float3_t* result);
void float4_reflect(const float4_t* v, const float4_t* normal, float4_t* result);

/// @brief calculates the projection of two vectors, storing into result
void float2_project(const float2_t* a, const float2_t* b, float2_t* result);
void float3_project(const float3_t* a, const float3_t* b, float3_t* result);
void float4_project(const float4_t* a, const float4_t* b, float4_t* result);

#endif // VECMATH_VEC_OP_INCLUDED