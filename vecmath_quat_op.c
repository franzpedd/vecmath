#include "vecmath_quat_op.h"
#include "vecmath_mat_op.h"
#include "vecmath_util.h"
#include <math.h>

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// identity
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API fquat fquat_identity()
{
    fquat result = { 0 };
    result.vector.x = 0;
    result.vector.y = 0;
    result.vector.z = 0;
    result.vector.w = 1;
    return result;
}

VECMATH_API dquat dquat_identity()
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

VECMATH_API fquat fquat_mul(const fquat* q1, const fquat* q2)
{
    fquat result = { 0 };
    result.vector.x = q1->vector.w * q2->vector.x + q1->vector.x * q2->vector.w + q1->vector.y * q2->vector.z - q1->vector.z * q2->vector.y;
    result.vector.y = q1->vector.w * q2->vector.y - q1->vector.x * q2->vector.z + q1->vector.y * q2->vector.w + q1->vector.z * q2->vector.x;
    result.vector.z = q1->vector.w * q2->vector.z + q1->vector.x * q2->vector.y - q1->vector.y * q2->vector.x + q1->vector.z * q2->vector.w;
    result.vector.w = q1->vector.w * q2->vector.w - q1->vector.x * q2->vector.x - q1->vector.y * q2->vector.y - q1->vector.z * q2->vector.z;
    return result;
}

VECMATH_API dquat dquat_mul(const dquat* q1, const dquat* q2)
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

VECMATH_API float fquat_length(const fquat* q)
{
    return sqrtf(q->vector.x * q->vector.x + q->vector.y * q->vector.y + q->vector.z * q->vector.z + q->vector.w * q->vector.w);
}

VECMATH_API double dquat_length(const dquat* q)
{
    return sqrt(q->vector.x * q->vector.x + q->vector.y * q->vector.y + q->vector.z * q->vector.z + q->vector.w * q->vector.w);
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// conjugate
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API fquat fquat_conjugate(const fquat* q)
{
    fquat result = { 0 };
    result.vector.x = -q->vector.x;
    result.vector.y = -q->vector.y;
    result.vector.z = -q->vector.z;
    result.vector.w = -q->vector.w;
    return result;
}

VECMATH_API dquat dquat_conjugate(const dquat* q)
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

VECMATH_API fquat fquat_normalize(const fquat *q)
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

VECMATH_API dquat dquat_normalize(const dquat* q)
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

VECMATH_API float fquat_dot(const fquat* q1, const fquat* q2)
{
    return q1->vector.x * q2->vector.x + q1->vector.y * q2->vector.y + q1->vector.z * q2->vector.z + q1->vector.w * q2->vector.w;
}

VECMATH_API double dquat_dot(const dquat* q1, const dquat* q2)
{
    return q1->vector.x * q2->vector.x + q1->vector.y * q2->vector.y + q1->vector.z * q2->vector.z + q1->vector.w * q2->vector.w;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// lerp
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API fquat fquat_lerp(const fquat *q1, const fquat *q2, float t)
{
    fquat result = { 0 };
    result.vector.x = q1->vector.x + t * (q2->vector.x - q1->vector.x);
    result.vector.y = q1->vector.y + t * (q2->vector.y - q1->vector.y);
    result.vector.z = q1->vector.z + t * (q2->vector.z - q1->vector.z);
    result.vector.w = q1->vector.w + t * (q2->vector.w - q1->vector.w);
    return result;
}

VECMATH_API dquat dquat_lerp(const dquat *q1, const dquat *q2, double t)
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

VECMATH_API fquat fquat_slerp(const fquat* q1, const fquat* q2, float t)
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

VECMATH_API dquat dquat_slerp(const dquat* q1, const dquat* q2, double t)
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

VECMATH_API fquat fquat_from_euler(const float3* rad)
{
    float cy = cosf(rad->xyz.y * 0.5f);
    float sy = sinf(rad->xyz.y * 0.5f);
    float cp = cosf(rad->xyz.x * 0.5f);
    float sp = sinf(rad->xyz.x * 0.5f);
    float cr = cosf(rad->xyz.z * 0.5f);
    float sr = sinf(rad->xyz.z * 0.5f);

    fquat q = { 0 };
    q.vector.w = cr * cp * cy + sr * sp * sy;
    q.vector.x = sr * cp * cy - cr * sp * sy;
    q.vector.y = cr * sp * cy + sr * cp * sy;
    q.vector.z = cr * cp * sy - sr * sp * cy;
    return q;
}

VECMATH_API dquat dquat_from_euler(const double3* d)
{
    double cy = cos(d->xyz.y * 0.5);
    double sy = sin(d->xyz.y * 0.5);
    double cp = cos(d->xyz.x * 0.5);
    double sp = sin(d->xyz.x * 0.5);
    double cr = cos(d->xyz.z * 0.5);
    double sr = sin(d->xyz.z * 0.5);

    dquat q = { 0 };
    q.vector.w = cr * cp * cy + sr * sp * sy;
    q.vector.x = sr * cp * cy - cr * sp * sy;
    q.vector.y = cr * sp * cy + sr * cp * sy;
    q.vector.z = cr * cp * sy - sr * sp * cy;
    return q;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// to_euler
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API float3 fquat_to_euler(const fquat* q)
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
        angles.xyz.y = copysignf((float)VECMATH_EPSILON_PI / 2.0f, sinp); // use 90 degrees if out of range
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

VECMATH_API double3 dquat_to_euler(const dquat *q)
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

VECMATH_API fmat4 fquat_to_fmat4_rowmajor(const fquat *q)
{
    fquat q_norm = fquat_normalize(q);
    float x = q_norm.vector.x, y = q_norm.vector.y, z = q_norm.vector.z, w = q_norm.vector.w;
    
    fmat4 result = { 0 };
    result.matrix.m00 = 1.0f - 2.0f * y * y - 2.0f * z * z;
    result.matrix.m01 = 2.0f * x * y - 2.0f * w * z;
    result.matrix.m02 = 2.0f * x * z + 2.0f * w * y;
    
    result.matrix.m10 = 2.0f * x * y + 2.0f * w * z;
    result.matrix.m11 = 1.0f - 2.0f*x*x - 2.0f * z * z;
    result.matrix.m12 = 2.0f * y * z - 2.0f * w * x;
    
    result.matrix.m20 = 2.0f * x * z - 2.0f * w * y;
    result.matrix.m21 = 2.0f * y * z + 2.0f * w * x;
    result.matrix.m22 = 1.0f - 2.0f * x * x - 2.0f * y * y;
    
    result.matrix.m33 = 1.0f;
    return result;
}

VECMATH_API fmat4 fquat_to_fmat4_colmajor(const fquat *q)
{
    fquat q_norm = fquat_normalize(q);
    float x = q_norm.vector.x, y = q_norm.vector.y, z = q_norm.vector.z, w = q_norm.vector.w;

    fmat4 result = { 0 };
    result.matrix.m00 = 1.0f - 2.0f * y * y - 2.0f * z * z;
    result.matrix.m10 = 2.0f * x * y + 2.0f * w * z;
    result.matrix.m20 = 2.0f * x * z - 2.0f * w * y;

    result.matrix.m01 = 2.0f * x * y - 2.0f * w * z;
    result.matrix.m11 = 1.0f - 2.0f * x * x - 2.0f * z * z;
    result.matrix.m21 = 2.0f * y * z + 2.0f * w * x;

    result.matrix.m02 = 2.0f * x * z + 2.0f * w * y;
    result.matrix.m12 = 2.0f * y * z - 2.0f * w * x;
    result.matrix.m22 = 1.0f - 2.0f * x * x - 2.0f * y * y;

    result.matrix.m33 = 1.0f;
    return result;
}

VECMATH_API dmat4 dquat_to_dmat4_rowmajor(const dquat *q)
{
    dquat q_norm = dquat_normalize(q);
    double x = q_norm.vector.x, y = q_norm.vector.y, z = q_norm.vector.z, w = q_norm.vector.w;
    
    dmat4 result = { 0 };
    result.matrix.m00 = 1.0 - 2.0 * y * y - 2.0 * z * z;
    result.matrix.m01 = 2.0 * x * y - 2.0 * w * z;
    result.matrix.m02 = 2.0 * x * z + 2.0 * w * y;

    result.matrix.m10 = 2.0 * x * y + 2.0 * w * z;
    result.matrix.m11 = 1.0 - 2.0 * x * x - 2.0 * z * z;
    result.matrix.m12 = 2.0 * y * z - 2.0 * w * x;

    result.matrix.m20 = 2.0 * x * z - 2.0 * w * y;
    result.matrix.m21 = 2.0 * y * z + 2.0 * w * x;
    result.matrix.m22 = 1.0 - 2.0 * x * x - 2.0 * y * y;

    result.matrix.m33 = 1.0;
    return result;
}

VECMATH_API dmat4 dquat_to_dmat4_colmajor(const dquat *q)
{
    dquat q_norm = dquat_normalize(q);
    double x = q_norm.vector.x, y = q_norm.vector.y, z = q_norm.vector.z, w = q_norm.vector.w;
    
    dmat4 result = { 0 };
    result.matrix.m00 = 1.0 - 2.0 * y * y - 2.0 * z * z;
    result.matrix.m10 = 2.0 * x * y + 2.0 * w * z;
    result.matrix.m20 = 2.0 * x * z - 2.0 * w * y;

    result.matrix.m01 = 2.0 * x * y - 2.0 * w * z;
    result.matrix.m11 = 1.0 - 2.0 * x * x - 2.0 * z * z;
    result.matrix.m21 = 2.0 * y * z + 2.0 * w * x;

    result.matrix.m02 = 2.0 * x * z + 2.0 * w * y;
    result.matrix.m12 = 2.0 * y * z - 2.0 * w * x;
    result.matrix.m22 = 1.0 - 2.0 * x * x - 2.0 * y * y;
    
    result.matrix.m33 = 1.0;
    return result;
}
