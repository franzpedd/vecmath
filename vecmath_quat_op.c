#include "vecmath_quat_op.h"

#include "vecmath_macros.h"
#include <math.h>

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// identity
/////////////////////////////////////////////////////////////////////////////////////

fquat fquat_identity()
{
    fquat result = { 0 };
    result.vector.x = 0;
    result.vector.y = 0;
    result.vector.z = 0;
    result.vector.w = 1;
    return result;
}

dquat dquat_identity()
{
    dquat result = { 0 };
    result.vector.x = 0;
    result.vector.y = 0;
    result.vector.z = 0;
    result.vector.w = 1;
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// mul
/////////////////////////////////////////////////////////////////////////////////////

fquat fquat_mul(const fquat* q1, const fquat* q2)
{
    fquat result = { 0 };
    result.vector.x = q1->vector.w * q2->vector.x + q1->vector.x * q2->vector.w + q1->vector.y * q2->vector.z - q1->vector.z * q2->vector.y;
    result.vector.y = q1->vector.w * q2->vector.y - q1->vector.x * q2->vector.z + q1->vector.y * q2->vector.w + q1->vector.z * q2->vector.x;
    result.vector.z = q1->vector.w * q2->vector.z + q1->vector.x * q2->vector.y - q1->vector.y * q2->vector.x + q1->vector.z * q2->vector.w;
    result.vector.w = q1->vector.w * q2->vector.w - q1->vector.x * q2->vector.x - q1->vector.y * q2->vector.y - q1->vector.z * q2->vector.z;
    return result;
}

dquat dquat_mul(const dquat* q1, const dquat* q2)
{
    dquat result = { 0 };
    result.vector.x = q1->vector.w * q2->vector.x + q1->vector.x * q2->vector.w + q1->vector.y * q2->vector.z - q1->vector.z * q2->vector.y;
    result.vector.y = q1->vector.w * q2->vector.y - q1->vector.x * q2->vector.z + q1->vector.y * q2->vector.w + q1->vector.z * q2->vector.x;
    result.vector.z = q1->vector.w * q2->vector.z + q1->vector.x * q2->vector.y - q1->vector.y * q2->vector.x + q1->vector.z * q2->vector.w;
    result.vector.w = q1->vector.w * q2->vector.w - q1->vector.x * q2->vector.x - q1->vector.y * q2->vector.y - q1->vector.z * q2->vector.z;
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// length
/////////////////////////////////////////////////////////////////////////////////////

float fquat_length(const fquat* q)
{
    return sqrtf(q->vector.x * q->vector.x + q->vector.y * q->vector.y + q->vector.z * q->vector.z + q->vector.w * q->vector.w);
}

double dquat_length(const dquat* q)
{
    return sqrt(q->vector.x * q->vector.x + q->vector.y * q->vector.y + q->vector.z * q->vector.z + q->vector.w * q->vector.w);
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// conjugate
/////////////////////////////////////////////////////////////////////////////////////

fquat fquat_conjugate(const fquat* q)
{
    fquat result = { 0 };
    result.vector.x = -q->vector.x;
    result.vector.y = -q->vector.y;
    result.vector.z = -q->vector.z;
    result.vector.w = -q->vector.w;
    return result;
}

dquat dquat_conjugate(const dquat* q)
{
    dquat result = { 0 };
    result.vector.x = -q->vector.x;
    result.vector.y = -q->vector.y;
    result.vector.z = -q->vector.z;
    result.vector.w = -q->vector.w;
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// normalize
/////////////////////////////////////////////////////////////////////////////////////

fquat fquat_normalize(const fquat *q)
{
    float len = sqrtf(q->vector.x * q->vector.x + q->vector.y * q->vector.y + q->vector.z * q->vector.z + q->vector.w * q->vector.w);
    if (len > VECMATH_EPSILON_FZERO) {
        fquat result = { 0 };
        result.vector.x = q->vector.x / len;
        result.vector.y = q->vector.y / len;
        result.vector.z = q->vector.z / len;
        result.vector.w = q->vector.w / len;
        return result;
    }
    return fquat_identity();
}

dquat dquat_normalize(const dquat* q)
{
    double len = sqrt(q->vector.x * q->vector.x + q->vector.y * q->vector.y + q->vector.z * q->vector.z + q->vector.w * q->vector.w);
    if (len > VECMATH_EPSILON_DZERO) {
        dquat result = { 0 };
        result.vector.x = q->vector.x / len;
        result.vector.y = q->vector.y / len;
        result.vector.z = q->vector.z / len;
        result.vector.w = q->vector.w / len;
        return result;
    }
    return dquat_identity();
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// dot
/////////////////////////////////////////////////////////////////////////////////////

float fquat_dot(const fquat* q1, const fquat* q2)
{
    return q1->vector.x * q2->vector.x + q1->vector.y * q2->vector.y + q1->vector.z * q2->vector.z + q1->vector.w * q2->vector.w;
}

double dquat_dot(const dquat* q1, const dquat* q2)
{
    return q1->vector.x * q2->vector.x + q1->vector.y * q2->vector.y + q1->vector.z * q2->vector.z + q1->vector.w * q2->vector.w;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// lerp
/////////////////////////////////////////////////////////////////////////////////////

fquat fquat_lerp(const fquat *q1, const fquat *q2, float t)
{
    fquat result = { 0 };
    result.vector.x = q1->vector.x + t * (q2->vector.x - q1->vector.x);
    result.vector.y = q1->vector.y + t * (q2->vector.y - q1->vector.y);
    result.vector.z = q1->vector.z + t * (q2->vector.z - q1->vector.z);
    result.vector.w = q1->vector.w + t * (q2->vector.w - q1->vector.w);
    return result;
}

dquat dquat_lerp(const dquat *q1, const dquat *q2, double t)
{
    dquat result = { 0 };
    result.vector.x = q1->vector.x + t * (q2->vector.x - q1->vector.x);
    result.vector.y = q1->vector.y + t * (q2->vector.y - q1->vector.y);
    result.vector.z = q1->vector.z + t * (q2->vector.z - q1->vector.z);
    result.vector.w = q1->vector.w + t * (q2->vector.w - q1->vector.w);
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// slerp
/////////////////////////////////////////////////////////////////////////////////////

fquat fquat_slerp(const fquat* q1, const fquat* q2, float t)
{
    float cos_theta = fquat_dot(q1, q2);
    
    if (fabsf(cos_theta) > 0.9999f) {
        return *q1;
    }
    
    fquat q2_adj = *q2;
    if (cos_theta < 0.0f) {
        q2_adj.vector.x = -q2->vector.x;
        q2_adj.vector.y = -q2->vector.y;
        q2_adj.vector.z = -q2->vector.z;
        q2_adj.vector.w = -q2->vector.w;
        cos_theta = -cos_theta;
    }
    
    float theta = acosf(cos_theta);
    float sin_theta = sqrtf(1.0f - cos_theta*cos_theta);
    
    if (fabsf(sin_theta) < VECMATH_EPSILON_FZERO) {
        return *q1;
    }
    
    float wa = sinf((1.0f - t) * theta) / sin_theta;
    float wb = sinf(t * theta) / sin_theta;
    
    fquat result = { 0 };
    result.vector.x = wa * q1->vector.x + wb * q2_adj.vector.x;
    result.vector.y = wa * q1->vector.y + wb * q2_adj.vector.y;
    result.vector.z = wa * q1->vector.z + wb * q2_adj.vector.z;
    result.vector.w = wa * q1->vector.w + wb * q2_adj.vector.w;
    return result;
}

dquat dquat_slerp(const dquat* q1, const dquat* q2, double t)
{
    double cos_theta = dquat_dot(q1, q2);
    
    if (fabs(cos_theta) > 0.9999) {
        return *q1;
    }
    
    dquat q2_adj = *q2;
    if (cos_theta < 0.0f) {
        q2_adj.vector.x = -q2->vector.x;
        q2_adj.vector.y = -q2->vector.y;
        q2_adj.vector.z = -q2->vector.z;
        q2_adj.vector.w = -q2->vector.w;
        cos_theta = -cos_theta;
    }
    
    double theta = acos(cos_theta);
    double sin_theta = sqrt(1.0 - cos_theta*cos_theta);
    
    if (fabs(sin_theta) < VECMATH_EPSILON_DZERO) {
        return *q1;
    }
    
    double wa = sin((1.0 - t) * theta) / sin_theta;
    double wb = sin(t * theta) / sin_theta;
    
    dquat result = { 0 };
    result.vector.x = wa * q1->vector.x + wb * q2_adj.vector.x;
    result.vector.y = wa * q1->vector.y + wb * q2_adj.vector.y;
    result.vector.z = wa * q1->vector.z + wb * q2_adj.vector.z;
    result.vector.w = wa * q1->vector.w + wb * q2_adj.vector.w;
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// from_euler
/////////////////////////////////////////////////////////////////////////////////////

fquat fquat_from_euler(const float3* axis, float angle_rad)
{
    float half_angle = angle_rad * 0.5f;
    float sin_half = sinf(half_angle);
    
    float len = sqrtf(axis->xyz.x * axis->xyz.x + axis->xyz.y * axis->xyz.y + axis->xyz.z * axis->xyz.z);
    if (len > VECMATH_EPSILON_FZERO) {
        float inv_len = 1.0f / len;
        fquat result = { 0 };
        result.vector.x = axis->xyz.x * inv_len * sin_half;
        result.vector.y = axis->xyz.y * inv_len * sin_half;
        result.vector.z = axis->xyz.z * inv_len * sin_half;
        result.vector.w = cosf(half_angle);
        return result;
    }
    return fquat_identity();
}

dquat dquat_from_euler(const double3 *axis, double angle_rad)
{
    double half_angle = angle_rad * 0.5f;
    double sin_half = sin(half_angle);
    
    double len = sqrt(axis->xyz.x * axis->xyz.x + axis->xyz.y * axis->xyz.y + axis->xyz.z * axis->xyz.z);
    if (len > VECMATH_EPSILON_DZERO) {
        double inv_len = 1.0 / len;
        dquat result = { 0 };
        result.vector.x = axis->xyz.x * inv_len * sin_half;
        result.vector.y = axis->xyz.y * inv_len * sin_half;
        result.vector.z = axis->xyz.z * inv_len * sin_half;
        result.vector.w = cos(half_angle);
        return result;
    }
    return dquat_identity();
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// to_euler
/////////////////////////////////////////////////////////////////////////////////////

float3 fquat_to_euler(const fquat* q)
{
    float3 angles = { 0 };
    fquat q_norm = fquat_normalize(q);
    float x = q_norm.vector.x, y = q_norm.vector.y, z = q_norm.vector.z, w = q_norm.vector.w;
    
    // roll (x-axis rotation)
    float sinr_cosp = 2.0f * (w * x + y * z);
    float cosr_cosp = 1.0f - 2.0f * (x * x + y * y);
    angles.xyz.x = atan2f(sinr_cosp, cosr_cosp);
    
    // pitch (y-axis rotation)
    float sinp = 2.0f * (w * y - z * x);
    if (fabsf(sinp) >= 1.0f) {
        angles.xyz.y = copysignf(VECMATH_EPSILON_PI / 2.0f, sinp); // use 90 degrees if out of range
    }

    else {
        angles.xyz.y = asinf(sinp);
    }
    
    // yaw (z-axis rotation)
    float siny_cosp = 2.0f * (w * z + x * y);
    float cosy_cosp = 1.0f - 2.0f * (y * y + z * z);
    angles.xyz.z = atan2f(siny_cosp, cosy_cosp);
    return angles;
}

double3 dquat_to_euler(const dquat *q)
{
    double3 angles = { 0 };
    dquat q_norm = dquat_normalize(q);
    double x = q_norm.vector.x, y = q_norm.vector.y, z = q_norm.vector.z, w = q_norm.vector.w;
    
    // roll (x-axis rotation)
    double sinr_cosp = 2.0 * (w * x + y * z);
    double cosr_cosp = 1.0 - 2.0 * (x * x + y * y);
    angles.xyz.x = atan2(sinr_cosp, cosr_cosp);
    
    // pitch (y-axis rotation)
    double sinp = 2.0 * (w * y - z * x);
    if (fabs(sinp) >= 1.0) {
        angles.xyz.y = copysign(VECMATH_EPSILON_PI / 2.0, sinp); // use 90 degrees if out of range
    }

    else {
        angles.xyz.y = asin(sinp);
    }
    
    // yaw (z-axis rotation)
    double siny_cosp = 2.0 * (w * z + x * y);
    double cosy_cosp = 1.0 - 2.0 * (y * y + z * z);
    angles.xyz.z = atan2(siny_cosp, cosy_cosp);
    return angles;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// to_matrix
/////////////////////////////////////////////////////////////////////////////////////

fmat4 fquat_to_fmat4_rowmajor(const fquat *q)
{
    fquat q_norm = fquat_normalize(q);
    float x = q_norm.vector.x, y = q_norm.vector.y, z = q_norm.vector.z, w = q_norm.vector.w;
    
    fmat4 result = { 0 };
    result.matrix.m00 = 1.0f - 2.0f*y*y - 2.0f*z*z;
    result.matrix.m01 = 2.0f*x*y - 2.0f*w*z;
    result.matrix.m02 = 2.0f*x*z + 2.0f*w*y;
    
    result.matrix.m10 = 2.0f*x*y + 2.0f*w*z;
    result.matrix.m11 = 1.0f - 2.0f*x*x - 2.0f*z*z;
    result.matrix.m12 = 2.0f*y*z - 2.0f*w*x;
    
    result.matrix.m20 = 2.0f*x*z - 2.0f*w*y;
    result.matrix.m21 = 2.0f*y*z + 2.0f*w*x;
    result.matrix.m22 = 1.0f - 2.0f*x*x - 2.0f*y*y;
    
    result.matrix.m33 = 1.0f;
    return result;
}

fmat4 fquat_to_fmat4_colmajor(const fquat *q)
{
    fquat q_norm = fquat_normalize(q);
    float x = q_norm.vector.x, y = q_norm.vector.y, z = q_norm.vector.z, w = q_norm.vector.w;
    
    fmat4 result = { 0 };
    result.matrix.m00 = 1.0f - 2.0f*y*y - 2.0f*z*z;
    result.matrix.m10 = 2.0f*x*y + 2.0f*w*z;
    result.matrix.m20 = 2.0f*x*z - 2.0f*w*y;
    
    result.matrix.m01 = 2.0f*x*y - 2.0f*w*z;
    result.matrix.m11 = 1.0f - 2.0f*x*x - 2.0f*z*z;
    result.matrix.m21 = 2.0f*y*z + 2.0f*w*x;
    
    result.matrix.m02 = 2.0f*x*z + 2.0f*w*y;
    result.matrix.m12 = 2.0f*y*z - 2.0f*w*x;
    result.matrix.m22 = 1.0f - 2.0f*x*x - 2.0f*y*y;
    
    result.matrix.m33 = 1.0f;
    return result;
}

dmat4 dquat_to_dmat4_rowmajor(const dquat *q)
{
    dquat q_norm = dquat_normalize(q);
    double x = q_norm.vector.x, y = q_norm.vector.y, z = q_norm.vector.z, w = q_norm.vector.w;
    
    dmat4 result = { 0 };
    result.matrix.m00 = 1.0 - 2.0*y*y - 2.0*z*z;
    result.matrix.m01 = 2.0*x*y - 2.0*w*z;
    result.matrix.m02 = 2.0*x*z + 2.0*w*y;
    
    result.matrix.m10 = 2.0*x*y + 2.0*w*z;
    result.matrix.m11 = 1.0 - 2.0*x*x - 2.0*z*z;
    result.matrix.m12 = 2.0*y*z - 2.0*w*x;
    
    result.matrix.m20 = 2.0*x*z - 2.0*w*y;
    result.matrix.m21 = 2.0*y*z + 2.0*w*x;
    result.matrix.m22 = 1.0 - 2.0*x*x - 2.0*y*y;
    
    result.matrix.m33 = 1.0;
    return result;
}

dmat4 dquat_to_dmat4_colmajor(const dquat *q)
{
    dquat q_norm = dquat_normalize(q);
    double x = q_norm.vector.x, y = q_norm.vector.y, z = q_norm.vector.z, w = q_norm.vector.w;
    
    dmat4 result = { 0 };
    result.matrix.m00 = 1.0 - 2.0*y*y - 2.0*z*z;
    result.matrix.m10 = 2.0*x*y + 2.0*w*z;
    result.matrix.m20 = 2.0*x*z - 2.0*w*y;
    
    result.matrix.m01 = 2.0*x*y - 2.0*w*z;
    result.matrix.m11 = 1.0 - 2.0*x*x - 2.0*z*z;
    result.matrix.m21 = 2.0*y*z + 2.0*w*x;
    
    result.matrix.m02 = 2.0*x*z + 2.0*w*y;
    result.matrix.m12 = 2.0*y*z - 2.0*w*x;
    result.matrix.m22 = 1.0 - 2.0*x*x - 2.0*y*y;
    
    result.matrix.m33 = 1.0;
    return result;
}
