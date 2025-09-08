#include "vecmath_vec_op.h"

#include <math.h>

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// length
/////////////////////////////////////////////////////////////////////////////////////

float float2_length(const float2* v)
{
    return sqrtf((v->xy.x * v->xy.x) + (v->xy.y * v->xy.y));
}

float float3_length(const float3* v)
{
    return sqrtf((v->xyz.x * v->xyz.x) + (v->xyz.y * v->xyz.y) + (v->xyz.z * v->xyz.z));
}

float float4_length(const float4* v)
{
    return sqrtf((v->xyzw.x * v->xyzw.x) + (v->xyzw.y * v->xyzw.y) + (v->xyzw.z * v->xyzw.z) + (v->xyzw.w * v->xyzw.w));
}

double double2_length(const double2* v)
{
    return sqrt((v->xy.x * v->xy.x) + (v->xy.y * v->xy.y));
}

double double3_length(const double3* v)
{
    return sqrt((v->xyz.x * v->xyz.x) + (v->xyz.y * v->xyz.y) + (v->xyz.z * v->xyz.z));
}

double double4_length(const double4* v)
{
    return sqrt((v->xyzw.x * v->xyzw.x) + (v->xyzw.y * v->xyzw.y) + (v->xyzw.z * v->xyzw.z) + (v->xyzw.w * v->xyzw.w));
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// length_sqrt 
/////////////////////////////////////////////////////////////////////////////////////

float float2_length_sqrt(const float2* v)
{
    return (v->xy.x * v->xy.x) + (v->xy.y * v->xy.y);
}

float float3_length_sqrt(const float3* v)
{
    return (v->xyz.x * v->xyz.x) + (v->xyz.y * v->xyz.y) + (v->xyz.z * v->xyz.z);
}

float float4_length_sqrt(const float4* v)
{
    return (v->xyzw.x * v->xyzw.x) + (v->xyzw.y * v->xyzw.y) + (v->xyzw.z * v->xyzw.z) + (v->xyzw.w * v->xyzw.w);
}

double double2_length_sqrt(const double2 *v)
{
    return (v->xy.x * v->xy.x) + (v->xy.y * v->xy.y);
}

double double3_length_sqrt(const double3 *v)
{
    return (v->xyz.x * v->xyz.x) + (v->xyz.y * v->xyz.y) + (v->xyz.z * v->xyz.z);
}

double double4_length_sqrt(const double4* v)
{
    return (v->xyzw.x * v->xyzw.x) + (v->xyzw.y * v->xyzw.y) + (v->xyzw.z * v->xyzw.z) + (v->xyzw.w * v->xyzw.w);
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// normalize 
/////////////////////////////////////////////////////////////////////////////////////

void float2_normalize(const float2* v, float2* result)
{
    float len = float2_length(v);
    if (len > 0.0f) {
        result->xy.x = v->xy.x / len;
        result->xy.y = v->xy.y / len;
        return;
    }
    result->xy.x = 0.0f;
    result->xy.y = 0.0f;
}

void float3_normalize(const float3* v, float3* result)
{
    float len = float3_length(v);
    if (len > 0.0f) {
        result->xyz.x = v->xyz.x / len;
        result->xyz.y = v->xyz.y / len;
        result->xyz.z = v->xyz.z / len;
        return;
    }
    result->xyz.x = 0.0f;
    result->xyz.y = 0.0f;
    result->xyz.z = 0.0f;
}

void float4_normalize(const float4* v, float4* result)
{
    float len = float4_length(v);
    if (len > 0.0f) {
        result->xyzw.x = v->xyzw.x / len;
        result->xyzw.y = v->xyzw.y / len;
        result->xyzw.z = v->xyzw.z / len;
        result->xyzw.w = v->xyzw.w / len;
        return;
    }
    result->xyzw.x = 0.0f;
    result->xyzw.y = 0.0f;
    result->xyzw.z = 0.0f;
    result->xyzw.w = 0.0f;
}

void double2_normalize(const double2* v, double2* result)
{
    double len = double2_length(v);
    if (len > 0.0f) {
        result->xy.x = v->xy.x / len;
        result->xy.y = v->xy.y / len;
        return;
    }
    result->xy.x = 0.0f;
    result->xy.y = 0.0f;
}

void double3_normalize(const double3 *v, double3* result)
{
    double len = double3_length(v);
    if (len > 0.0f) {
        result->xyz.x = v->xyz.x / len;
        result->xyz.y = v->xyz.y / len;
        result->xyz.z = v->xyz.z / len;
        return;
    }
    result->xyz.x = 0.0f;
    result->xyz.y = 0.0f;
    result->xyz.z = 0.0f;
}

void double4_normalize(const double4* v, double4* result)
{
    double len = double4_length(v);
    if (len > 0.0f) {
        result->xyzw.x = v->xyzw.x / len;
        result->xyzw.y = v->xyzw.y / len;
        result->xyzw.z = v->xyzw.z / len;
        result->xyzw.w = v->xyzw.w / len;
        return;
    }
    result->xyzw.x = 0.0f;
    result->xyzw.y = 0.0f;
    result->xyzw.z = 0.0f;
    result->xyzw.w = 0.0f;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// dot 
/////////////////////////////////////////////////////////////////////////////////////

float float2_dot(const float2 *a, const float2 *b)
{
    return (a->xy.x * b->xy.x) + (a->xy.y * b->xy.y);
}

float float3_dot(const float3* a, const float3* b)
{
    return (a->xyz.x * b->xyz.x) + (a->xyz.y * b->xyz.y) + (a->xyz.z * b->xyz.z);
}

float float4_dot(const float4* a, const float4* b)
{
    return (a->xyzw.x * b->xyzw.x) + (a->xyzw.y * b->xyzw.y) + (a->xyzw.z * b->xyzw.z) + (a->xyzw.w * b->xyzw.w);
}

double double2_dot(const double2* a, const double2* b)
{
    return (a->xy.x * b->xy.x) + (a->xy.y * b->xy.y);
}

double double3_dot(const double3* a, const double3* b)
{
    return (a->xyz.x * b->xyz.x) + (a->xyz.y * b->xyz.y) + (a->xyz.z * b->xyz.z);
}

double double4_dot(const double4* a, const double4* b)
{
    return (a->xyzw.x * b->xyzw.x) + (a->xyzw.y * b->xyzw.y) + (a->xyzw.z * b->xyzw.z) + (a->xyzw.w * b->xyzw.w);
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// cross 
/////////////////////////////////////////////////////////////////////////////////////

void float2_cross(const float2* a, const float2* b, float* result)
{
    *result = (a->xy.x * b->xy.y) - (a->xy.y * b->xy.x);
}

void float3_cross(const float3* a, const float3* b, float3* result)
{
    result->xyz.x = a->xyz.y * b->xyz.z - a->xyz.z * b->xyz.y;
    result->xyz.y = a->xyz.z * b->xyz.x - a->xyz.x * b->xyz.z;
    result->xyz.z = a->xyz.x * b->xyz.y - a->xyz.y * b->xyz.x;
}

void double2_cross(const double2* a, const double2* b, double* result)
{
    *result = (a->xy.x * b->xy.y) - (a->xy.y * b->xy.x);
}

void double3_cross(const double3* a, const double3* b, double3* result)
{
    result->xyz.x = a->xyz.y * b->xyz.z - a->xyz.z * b->xyz.y;
    result->xyz.y = a->xyz.z * b->xyz.x - a->xyz.x * b->xyz.z;
    result->xyz.z = a->xyz.x * b->xyz.y - a->xyz.y * b->xyz.x;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// scale 
/////////////////////////////////////////////////////////////////////////////////////

void float2_scale(const float2* v, float scalar, float2* result)
{
    result->xy.x = v->xy.x * scalar;
    result->xy.y = v->xy.y * scalar;
}

void float3_scale(const float3* v, float scalar, float3* result)
{
    result->xyz.x = v->xyz.x * scalar;
    result->xyz.y = v->xyz.y * scalar;
    result->xyz.z = v->xyz.z * scalar;
}

void float4_scale(const float4* v, float scalar, float4* result)
{
    result->xyzw.x = v->xyzw.x * scalar;
    result->xyzw.y = v->xyzw.y * scalar;
    result->xyzw.z = v->xyzw.z * scalar;
    result->xyzw.w = v->xyzw.w * scalar;
}

void double2_scale(const double2* v, double scalar, double2* result)
{
    result->xy.x = v->xy.x * scalar;
    result->xy.y = v->xy.y * scalar;
}

void double3_scale(const double3* v, double scalar, double3* result)
{
    result->xyz.x = v->xyz.x * scalar;
    result->xyz.y = v->xyz.y * scalar;
    result->xyz.z = v->xyz.z * scalar;
}

void double4_scale(const double4* v, double scalar, double4* result)
{
    result->xyzw.x = v->xyzw.x * scalar;
    result->xyzw.y = v->xyzw.y * scalar;
    result->xyzw.z = v->xyzw.z * scalar;
    result->xyzw.w = v->xyzw.w * scalar;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// lerp 
/////////////////////////////////////////////////////////////////////////////////////

void float2_lerp(const float2* a, const float2* b, float t, float2* result)
{
    result->xy.x = a->xy.x + (b->xy.x - a->xy.x) * t;
    result->xy.y = a->xy.y + (b->xy.y - a->xy.y) * t;
}

void float3_lerp(const float3* a, const float3* b, float t, float3* result)
{
    result->xyz.x = a->xyz.x + (b->xyz.x - a->xyz.x) * t;
    result->xyz.y = a->xyz.y + (b->xyz.y - a->xyz.y) * t;
    result->xyz.z = a->xyz.z + (b->xyz.z - a->xyz.z) * t;
}

void float4_lerp(const float4* a, const float4* b, float t, float4* result)
{
    result->xyzw.x = a->xyzw.x + (b->xyzw.x - a->xyzw.x) * t;
    result->xyzw.y = a->xyzw.y + (b->xyzw.y - a->xyzw.y) * t;
    result->xyzw.z = a->xyzw.z + (b->xyzw.z - a->xyzw.z) * t;
    result->xyzw.w = a->xyzw.w + (b->xyzw.w - a->xyzw.w) * t;
}

void double2_lerp(const double2* a, const double2* b, double t, double2* result)
{
    result->xy.x = a->xy.x + (b->xy.x - a->xy.x) * t;
    result->xy.y = a->xy.y + (b->xy.y - a->xy.y) * t;
}

void double3_lerp(const double3* a, const double3* b, double t, double3* result)
{
    result->xyz.x = a->xyz.x + (b->xyz.x - a->xyz.x) * t;
    result->xyz.y = a->xyz.y + (b->xyz.y - a->xyz.y) * t;
    result->xyz.z = a->xyz.z + (b->xyz.z - a->xyz.z) * t;
}

void double4_lerp(const double4* a, const double4* b, double t, double4* result)
{
    result->xyzw.x = a->xyzw.x + (b->xyzw.x - a->xyzw.x) * t;
    result->xyzw.y = a->xyzw.y + (b->xyzw.y - a->xyzw.y) * t;
    result->xyzw.z = a->xyzw.z + (b->xyzw.z - a->xyzw.z) * t;
    result->xyzw.w = a->xyzw.w + (b->xyzw.w - a->xyzw.w) * t;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// distance 
/////////////////////////////////////////////////////////////////////////////////////

float float2_distance(const float2* a, const float2* b)
{
    float dx = a->xy.x - b->xy.x;
    float dy = a->xy.y - b->xy.y;
    return sqrtf(dx * dx + dy * dy);
}

float float3_distance(const float3* a, const float3* b)
{
    float dx = a->xyz.x - b->xyz.x;
    float dy = a->xyz.y - b->xyz.y;
    float dz = a->xyz.z - b->xyz.z;
    return sqrtf(dx * dx + dy * dy + dz * dz);
}

float float4_distance(const float4* a, const float4* b)
{
    float dx = a->xyzw.x - b->xyzw.x;
    float dy = a->xyzw.y - b->xyzw.y;
    float dz = a->xyzw.z - b->xyzw.z;
    float dw = a->xyzw.w - b->xyzw.w;
    return sqrtf(dx * dx + dy * dy + dz * dz + dw * dw);
}

double double2_distance(const double2* a, const double2* b)
{
    double dx = a->xy.x - b->xy.x;
    double dy = a->xy.y - b->xy.y;
    return sqrt(dx * dx + dy * dy);
}

double double3_distance(const double3* a, const double3* b)
{
    double dx = a->xyz.x - b->xyz.x;
    double dy = a->xyz.y - b->xyz.y;
    double dz = a->xyz.z - b->xyz.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

double double4_distance(const double4* a, const double4* b)
{
    double dx = a->xyzw.x - b->xyzw.x;
    double dy = a->xyzw.y - b->xyzw.y;
    double dz = a->xyzw.z - b->xyzw.z;
    double dw = a->xyzw.w - b->xyzw.w;
    return sqrt(dx * dx + dy * dy + dz * dz + dw * dw);
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// distance_sqrt 
/////////////////////////////////////////////////////////////////////////////////////

float float2_distance_sqrt(const float2* a, const float2* b)
{
    float dx = a->xy.x - b->xy.x;
    float dy = a->xy.y - b->xy.y;
    return dx * dx + dy * dy;
}

float float3_distance_sqrt(const float3* a, const float3* b)
{
    float dx = a->xyz.x - b->xyz.x;
    float dy = a->xyz.y - b->xyz.y;
    float dz = a->xyz.z - b->xyz.z;
    return dx * dx + dy * dy + dz * dz;
}

float float4_distance_sqrt(const float4* a, const float4* b)
{
    float dx = a->xyzw.x - b->xyzw.x;
    float dy = a->xyzw.y - b->xyzw.y;
    float dz = a->xyzw.z - b->xyzw.z;
    float dw = a->xyzw.w - b->xyzw.w;
    return dx * dx + dy * dy + dz * dz + dw * dw;
}

double double2_distance_sqrt(const double2* a, const double2*b )
{
    double dx = a->xy.x - b->xy.x;
    double dy = a->xy.y - b->xy.y;
    return dx * dx + dy * dy;
}

double double3_distance_sqrt(const double3* a, const double3* b)
{
    double dx = a->xyz.x - b->xyz.x;
    double dy = a->xyz.y - b->xyz.y;
    double dz = a->xyz.z - b->xyz.z;
    return dx * dx + dy * dy + dz * dz;
}

double double4_distance_sqrt(const double4* a, const double4* b)
{
    double dx = a->xyzw.x - b->xyzw.x;
    double dy = a->xyzw.y - b->xyzw.y;
    double dz = a->xyzw.z - b->xyzw.z;
    double dw = a->xyzw.w - b->xyzw.w;
    return dx * dx + dy * dy + dz * dz + dw * dw;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// reflect
/////////////////////////////////////////////////////////////////////////////////////

void float2_reflect(const float2* v, const float2* normal, float2* result)
{
    double dot = float2_dot(v, normal);
    double normal_length_sq = float2_dot(normal, normal);
    
    if (normal_length_sq > 0.0f) {
        float factor = 2.0f * dot / normal_length_sq;
        result->xy.x = v->xy.x - factor * normal->xy.x;
        result->xy.y = v->xy.y - factor * normal->xy.y;
        return;
    }
    *result = *v;
}

void float3_reflect(const float3* v, const float3* normal, float3* result)
{
    float dot = float3_dot(v, normal);
    float normal_length_sq = float3_dot(normal, normal);
    
    if (normal_length_sq > 0.0f) {
        float factor = 2.0f * dot / normal_length_sq;
        result->xyz.x = v->xyz.x - factor * normal->xyz.x;
        result->xyz.y = v->xyz.y - factor * normal->xyz.y;
        result->xyz.z = v->xyz.z - factor * normal->xyz.z;
        return;
    }
    *result = *v;
}

void float4_reflect(const float4* v, const float4* normal, float4* result)
{
    float dot = float4_dot(v, normal);
    float normal_length_sq = float4_dot(normal, normal);
    
    if (normal_length_sq > 0.0f) {
        float factor = 2.0f * dot / normal_length_sq;
        result->xyzw.x = v->xyzw.x - factor * normal->xyzw.x;
        result->xyzw.y = v->xyzw.y - factor * normal->xyzw.y;
        result->xyzw.z = v->xyzw.z - factor * normal->xyzw.z;
        result->xyzw.w = v->xyzw.w - factor * normal->xyzw.w;
        return;
    }
    *result = *v;
}

void double2_reflect(const double2* v, const double2* normal, double2* result)
{
    double dot = double2_dot(v, normal);
    double normal_length_sq = double2_dot(normal, normal);
    
    if (normal_length_sq > 0.0f) {
        double factor = 2.0f * dot / normal_length_sq;
        result->xy.x = v->xy.x - factor * normal->xy.x;
        result->xy.y = v->xy.y - factor * normal->xy.y;
        return;
    }
    *result = *v;
}

void double3_reflect(const double3* v, const double3* normal, double3* result)
{
    double dot = double3_dot(v, normal);
    double normal_length_sq = double3_dot(normal, normal);
    
    if (normal_length_sq > 0.0f) {
        double factor = 2.0f * dot / normal_length_sq;
        result->xyz.x = v->xyz.x - factor * normal->xyz.x;
        result->xyz.y = v->xyz.y - factor * normal->xyz.y;
        result->xyz.z = v->xyz.z - factor * normal->xyz.z;
        return;
    }
    *result = *v;
}

void double4_reflect(const double4* v, const double4* normal, double4* result)
{
    double dot = double4_dot(v, normal);
    double normal_length_sq = double4_dot(normal, normal);
    
    if (normal_length_sq > 0.0f) {
        double factor = 2.0f * dot / normal_length_sq;
        result->xyzw.x = v->xyzw.x - factor * normal->xyzw.x;
        result->xyzw.y = v->xyzw.y - factor * normal->xyzw.y;
        result->xyzw.z = v->xyzw.z - factor * normal->xyzw.z;
        result->xyzw.w = v->xyzw.w - factor * normal->xyzw.w;
        return;
    }
    *result = *v;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// project
/////////////////////////////////////////////////////////////////////////////////////

void float2_project(const float2* a, const float2* b, float2* result)
{
    float dot_ab = float2_dot(a, b);
    float dot_bb = float2_dot(b, b);
    
    if (dot_bb > 0.0f) {
        float scale = dot_ab / dot_bb;
        result->xy.x = b->xy.x * scale;
        result->xy.y = b->xy.y * scale;
        return;
    }
    result->xy.x = 0.0f;
    result->xy.y = 0.0f;
}

void float3_project(const float3* a, const float3* b, float3* result)
{
    float dot_ab = float3_dot(a, b);
    float dot_bb = float3_dot(b, b);
    
    if (dot_bb > 0.0f) {
        float scale = dot_ab / dot_bb;
        result->xyz.x = b->xyz.x * scale;
        result->xyz.y = b->xyz.y * scale;
        result->xyz.z = b->xyz.z * scale;
        return;
    }
    result->xyz.x = 0.0f;
    result->xyz.y = 0.0f;
    result->xyz.z = 0.0f;
}

void float4_project(const float4 *a, const float4 *b, float4 *result)
{
    float dot_ab = float4_dot(a, b);
    float dot_bb = float4_dot(b, b);
    
    if (dot_bb > 0.0f) {
        float scale = dot_ab / dot_bb;
        result->xyzw.x = b->xyzw.x * scale;
        result->xyzw.y = b->xyzw.y * scale;
        result->xyzw.z = b->xyzw.z * scale;
        result->xyzw.w = b->xyzw.w * scale;
    }
    result->xyzw.x = 0.0f;
    result->xyzw.y = 0.0f;
    result->xyzw.z = 0.0f;
    result->xyzw.w = 0.0f;
}

void double2_project(const double2* a, const double2* b, double2* result)
{
    double dot_ab = double2_dot(a, b);
    double dot_bb = double2_dot(b, b);
    
    if (dot_bb > 0.0f) {
        double scale = dot_ab / dot_bb;
        result->xy.x = b->xy.x * scale;
        result->xy.y = b->xy.y * scale;
        return;
    }
    result->xy.x = 0.0f;
    result->xy.y = 0.0f;
}

void double3_project(const double3 *a, const double3 *b, double3 *result)
{
    double dot_ab = double3_dot(a, b);
    double dot_bb = double3_dot(b, b);
    
    if (dot_bb > 0.0f) {
        double scale = dot_ab / dot_bb;
        result->xyz.x = b->xyz.x * scale;
        result->xyz.y = b->xyz.y * scale;
        result->xyz.z = b->xyz.z * scale;
        return;
    }
    result->xyz.x = 0.0f;
    result->xyz.y = 0.0f;
    result->xyz.z = 0.0f;
}

void double4_project(const double4 *a, const double4 *b, double4 *result)
{
    double dot_ab = double4_dot(a, b);
    double dot_bb = double4_dot(b, b);
    
    if (dot_bb > 0.0f) {
        double scale = dot_ab / dot_bb;
        result->xyzw.x = b->xyzw.x * scale;
        result->xyzw.y = b->xyzw.y * scale;
        result->xyzw.z = b->xyzw.z * scale;
        result->xyzw.w = b->xyzw.w * scale;
    }
    result->xyzw.x = 0.0f;
    result->xyzw.y = 0.0f;
    result->xyzw.z = 0.0f;
    result->xyzw.w = 0.0f;
}
