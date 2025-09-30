#include "vecmath_basic_op.h"
#include <math.h>
#include <string.h>

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// add 
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API float2 float2_add(const float2* a, const float2* b)
{
    float2 result = { 0 };
    result.xy.x = a->xy.x + b->xy.x;
    result.xy.y = a->xy.y + b->xy.y;
    return result;
}

VECMATH_API float3 float3_add(const float3* a, const float3* b)
{
    float3 result = { 0 };
    result.xyz.x = a->xyz.x + b->xyz.x;
    result.xyz.y = a->xyz.y + b->xyz.y;
    result.xyz.z = a->xyz.z + b->xyz.z;
    return result;
}

VECMATH_API float4 float4_add(const float4* a, const float4* b)
{
    float4 result = { 0 };
    result.xyzw.x = a->xyzw.x + b->xyzw.x;
    result.xyzw.y = a->xyzw.y + b->xyzw.y;
    result.xyzw.z = a->xyzw.z + b->xyzw.z;
    result.xyzw.w = a->xyzw.w + b->xyzw.w;
    return result;
}

VECMATH_API fmat2 fmat2_add(const fmat2* a, const fmat2* b)
{
    fmat2 result = { 0 };
    result.matrix.m00 = a->matrix.m00 + b->matrix.m00;
    result.matrix.m01 = a->matrix.m01 + b->matrix.m01;

    result.matrix.m10 = a->matrix.m10 + b->matrix.m10;
    result.matrix.m11 = a->matrix.m11 + b->matrix.m11;
    return result;
}

VECMATH_API fmat3 fmat3_add(const fmat3* a, const fmat3* b)
{
    fmat3 result = { 0 };
    result.matrix.m00 = a->matrix.m00 + b->matrix.m00;
    result.matrix.m01 = a->matrix.m01 + b->matrix.m01;
    result.matrix.m02 = a->matrix.m02 + b->matrix.m02;

    result.matrix.m10 = a->matrix.m10 + b->matrix.m10;
    result.matrix.m11 = a->matrix.m11 + b->matrix.m11;
    result.matrix.m12 = a->matrix.m12 + b->matrix.m12;

    result.matrix.m20 = a->matrix.m20 + b->matrix.m20;
    result.matrix.m21 = a->matrix.m21 + b->matrix.m21;
    result.matrix.m22 = a->matrix.m22 + b->matrix.m22;
    return result;
}

VECMATH_API fmat4 fmat4_add(const fmat4* a, const fmat4* b)
{
    fmat4 result = { 0 };
    result.matrix.m00 = a->matrix.m00 + b->matrix.m00;
    result.matrix.m01 = a->matrix.m01 + b->matrix.m01;
    result.matrix.m02 = a->matrix.m02 + b->matrix.m02;
    result.matrix.m03 = a->matrix.m03 + b->matrix.m03;

    result.matrix.m10 = a->matrix.m10 + b->matrix.m10;
    result.matrix.m11 = a->matrix.m11 + b->matrix.m11;
    result.matrix.m12 = a->matrix.m12 + b->matrix.m12;
    result.matrix.m13 = a->matrix.m13 + b->matrix.m13;

    result.matrix.m20 = a->matrix.m20 + b->matrix.m20;
    result.matrix.m21 = a->matrix.m21 + b->matrix.m21;
    result.matrix.m22 = a->matrix.m22 + b->matrix.m22;
    result.matrix.m23 = a->matrix.m23 + b->matrix.m23;

    result.matrix.m30 = a->matrix.m30 + b->matrix.m30;
    result.matrix.m31 = a->matrix.m31 + b->matrix.m31;
    result.matrix.m32 = a->matrix.m32 + b->matrix.m32;
    result.matrix.m33 = a->matrix.m33 + b->matrix.m33;
    return result;
}

VECMATH_API double2 double2_add(const double2* a, const double2* b)
{
    double2 result = { 0 };
    result.xy.x = a->xy.x + b->xy.x;
    result.xy.y = a->xy.y + b->xy.y;
    return result;
}

VECMATH_API double3 double3_add(const double3* a, const double3* b)
{
    double3 result = { 0 };
    result.xyz.x = a->xyz.x + b->xyz.x;
    result.xyz.y = a->xyz.y + b->xyz.y;
    result.xyz.z = a->xyz.z + b->xyz.z;
    return result;
}

VECMATH_API double4 double4_add(const double4* a, const double4* b)
{
    double4 result = { 0 };
    result.xyzw.x = a->xyzw.x + b->xyzw.x;
    result.xyzw.y = a->xyzw.y + b->xyzw.y;
    result.xyzw.z = a->xyzw.z + b->xyzw.z;
    result.xyzw.w = a->xyzw.w + b->xyzw.w;
    return result;
}

VECMATH_API dmat2 dmat2_add(const dmat2* a, const dmat2* b)
{
    dmat2 result = { 0 };
    result.matrix.m00 = a->matrix.m00 + b->matrix.m00;
    result.matrix.m01 = a->matrix.m01 + b->matrix.m01;

    result.matrix.m10 = a->matrix.m10 + b->matrix.m10;
    result.matrix.m11 = a->matrix.m11 + b->matrix.m11;
    return result;
}

VECMATH_API dmat3 dmat3_add(const dmat3* a, const dmat3* b)
{
    dmat3 result = { 0 };
    result.matrix.m00 = a->matrix.m00 + b->matrix.m00;
    result.matrix.m01 = a->matrix.m01 + b->matrix.m01;
    result.matrix.m02 = a->matrix.m02 + b->matrix.m02;

    result.matrix.m10 = a->matrix.m10 + b->matrix.m10;
    result.matrix.m11 = a->matrix.m11 + b->matrix.m11;
    result.matrix.m12 = a->matrix.m12 + b->matrix.m12;

    result.matrix.m20 = a->matrix.m20 + b->matrix.m20;
    result.matrix.m21 = a->matrix.m21 + b->matrix.m21;
    result.matrix.m22 = a->matrix.m22 + b->matrix.m22;
    return result;
}

VECMATH_API dmat4 dmat4_add(const dmat4* a, const dmat4* b)
{
    dmat4 result = { 0 };
    result.matrix.m00 = a->matrix.m00 + b->matrix.m00;
    result.matrix.m01 = a->matrix.m01 + b->matrix.m01;
    result.matrix.m02 = a->matrix.m02 + b->matrix.m02;
    result.matrix.m03 = a->matrix.m03 + b->matrix.m03;

    result.matrix.m10 = a->matrix.m10 + b->matrix.m10;
    result.matrix.m11 = a->matrix.m11 + b->matrix.m11;
    result.matrix.m12 = a->matrix.m12 + b->matrix.m12;
    result.matrix.m13 = a->matrix.m13 + b->matrix.m13;

    result.matrix.m20 = a->matrix.m20 + b->matrix.m20;
    result.matrix.m21 = a->matrix.m21 + b->matrix.m21;
    result.matrix.m22 = a->matrix.m22 + b->matrix.m22;
    result.matrix.m23 = a->matrix.m23 + b->matrix.m23;

    result.matrix.m30 = a->matrix.m30 + b->matrix.m30;
    result.matrix.m31 = a->matrix.m31 + b->matrix.m31;
    result.matrix.m32 = a->matrix.m32 + b->matrix.m32;
    result.matrix.m33 = a->matrix.m33 + b->matrix.m33;
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// sub
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API float2 float2_sub(const float2* a, const float2* b)
{
    float2 result = { 0 };
    result.xy.x = a->xy.x - b->xy.x;
    result.xy.y = a->xy.y - b->xy.y;
    return result;
}

VECMATH_API float3 float3_sub(const float3* a, const float3* b)
{
    float3 result = { 0 };
    result.xyz.x = a->xyz.x - b->xyz.x;
    result.xyz.y = a->xyz.y - b->xyz.y;
    result.xyz.z = a->xyz.z - b->xyz.z;
    return result;
}

VECMATH_API float4 float4_sub(const float4* a, const float4* b)
{
    float4 result = { 0 };
    result.xyzw.x = a->xyzw.x - b->xyzw.x;
    result.xyzw.y = a->xyzw.y - b->xyzw.y;
    result.xyzw.z = a->xyzw.z - b->xyzw.z;
    result.xyzw.w = a->xyzw.w - b->xyzw.w;
    return result;
}

VECMATH_API fmat2 fmat2_sub(const fmat2* a, const fmat2* b)
{
    fmat2 result = { 0 };
    result.matrix.m00 = a->matrix.m00 - b->matrix.m00;
    result.matrix.m01 = a->matrix.m01 - b->matrix.m01;

    result.matrix.m10 = a->matrix.m10 - b->matrix.m10;
    result.matrix.m11 = a->matrix.m11 - b->matrix.m11;
    return result;
}

VECMATH_API fmat3 fmat3_sub(const fmat3* a, const fmat3* b)
{
    fmat3 result = { 0 };
    result.matrix.m00 = a->matrix.m00 - b->matrix.m00;
    result.matrix.m01 = a->matrix.m01 - b->matrix.m01;
    result.matrix.m02 = a->matrix.m02 - b->matrix.m02;

    result.matrix.m10 = a->matrix.m10 - b->matrix.m10;
    result.matrix.m11 = a->matrix.m11 - b->matrix.m11;
    result.matrix.m12 = a->matrix.m12 - b->matrix.m12;

    result.matrix.m20 = a->matrix.m20 - b->matrix.m20;
    result.matrix.m21 = a->matrix.m21 - b->matrix.m21;
    result.matrix.m22 = a->matrix.m22 - b->matrix.m22;
    return result;
}

VECMATH_API fmat4 fmat4_sub(const fmat4* a, const fmat4* b)
{
    fmat4 result = { 0 };
    result.matrix.m00 = a->matrix.m00 - b->matrix.m00;
    result.matrix.m01 = a->matrix.m01 - b->matrix.m01;
    result.matrix.m02 = a->matrix.m02 - b->matrix.m02;
    result.matrix.m03 = a->matrix.m03 - b->matrix.m03;

    result.matrix.m10 = a->matrix.m10 - b->matrix.m10;
    result.matrix.m11 = a->matrix.m11 - b->matrix.m11;
    result.matrix.m12 = a->matrix.m12 - b->matrix.m12;
    result.matrix.m13 = a->matrix.m13 - b->matrix.m13;

    result.matrix.m20 = a->matrix.m20 - b->matrix.m20;
    result.matrix.m21 = a->matrix.m21 - b->matrix.m21;
    result.matrix.m22 = a->matrix.m22 - b->matrix.m22;
    result.matrix.m23 = a->matrix.m23 - b->matrix.m23;

    result.matrix.m30 = a->matrix.m30 - b->matrix.m30;
    result.matrix.m31 = a->matrix.m31 - b->matrix.m31;
    result.matrix.m32 = a->matrix.m32 - b->matrix.m32;
    result.matrix.m33 = a->matrix.m33 - b->matrix.m33;
    return result;
}

VECMATH_API double2 double2_sub(const double2* a, const double2* b)
{
    double2 result = { 0 };
    result.xy.x = a->xy.x - b->xy.x;
    result.xy.y = a->xy.y - b->xy.y;
    return result;
}

VECMATH_API double3 double3_sub(const double3* a, const double3* b)
{
    double3 result = { 0 };
    result.xyz.x = a->xyz.x - b->xyz.x;
    result.xyz.y = a->xyz.y - b->xyz.y;
    result.xyz.z = a->xyz.z - b->xyz.z;
    return result;
}

VECMATH_API double4 double4_sub(const double4* a, const double4* b)
{
    double4 result = { 0 };
    result.xyzw.x = a->xyzw.x - b->xyzw.x;
    result.xyzw.y = a->xyzw.y - b->xyzw.y;
    result.xyzw.z = a->xyzw.z - b->xyzw.z;
    result.xyzw.w = a->xyzw.w - b->xyzw.w;
    return result;
}

VECMATH_API dmat2 dmat2_sub(const dmat2* a, const dmat2* b)
{
    dmat2 result = { 0 };
    result.matrix.m00 = a->matrix.m00 - b->matrix.m00;
    result.matrix.m01 = a->matrix.m01 - b->matrix.m01;

    result.matrix.m10 = a->matrix.m10 - b->matrix.m10;
    result.matrix.m11 = a->matrix.m11 - b->matrix.m11;
    return result;
}

VECMATH_API dmat3 dmat3_sub(const dmat3* a, const dmat3* b)
{
    dmat3 result = { 0 };
    result.matrix.m00 = a->matrix.m00 - b->matrix.m00;
    result.matrix.m01 = a->matrix.m01 - b->matrix.m01;
    result.matrix.m02 = a->matrix.m02 - b->matrix.m02;

    result.matrix.m10 = a->matrix.m10 - b->matrix.m10;
    result.matrix.m11 = a->matrix.m11 - b->matrix.m11;
    result.matrix.m12 = a->matrix.m12 - b->matrix.m12;

    result.matrix.m20 = a->matrix.m20 - b->matrix.m20;
    result.matrix.m21 = a->matrix.m21 - b->matrix.m21;
    result.matrix.m22 = a->matrix.m22 - b->matrix.m22;
    return result;
}

VECMATH_API dmat4 dmat4_sub(const dmat4* a, const dmat4* b)
{
    dmat4 result = { 0 };
    result.matrix.m00 = a->matrix.m00 - b->matrix.m00;
    result.matrix.m01 = a->matrix.m01 - b->matrix.m01;
    result.matrix.m02 = a->matrix.m02 - b->matrix.m02;
    result.matrix.m03 = a->matrix.m03 - b->matrix.m03;

    result.matrix.m10 = a->matrix.m10 - b->matrix.m10;
    result.matrix.m11 = a->matrix.m11 - b->matrix.m11;
    result.matrix.m12 = a->matrix.m12 - b->matrix.m12;
    result.matrix.m13 = a->matrix.m13 - b->matrix.m13;

    result.matrix.m20 = a->matrix.m20 - b->matrix.m20;
    result.matrix.m21 = a->matrix.m21 - b->matrix.m21;
    result.matrix.m22 = a->matrix.m22 - b->matrix.m22;
    result.matrix.m23 = a->matrix.m23 - b->matrix.m23;

    result.matrix.m30 = a->matrix.m30 - b->matrix.m30;
    result.matrix.m31 = a->matrix.m31 - b->matrix.m31;
    result.matrix.m32 = a->matrix.m32 - b->matrix.m32;
    result.matrix.m33 = a->matrix.m33 - b->matrix.m33;
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// mul
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API float2 float2_mul(const float2* a, const float2* b)
{
    float2 result = { 0 };
    result.xy.x = a->xy.x * b->xy.x;
    result.xy.y = a->xy.y * b->xy.y;
    return result;
}

VECMATH_API float3 float3_mul(const float3* a, const float3* b)
{
    float3 result = { 0 };
    result.xyz.x = a->xyz.x * b->xyz.x;
    result.xyz.y = a->xyz.y * b->xyz.y;
    result.xyz.z = a->xyz.z * b->xyz.z;
    return result;
}

VECMATH_API float4 float4_mul(const float4* a, const float4* b)
{
    float4 result = { 0 };
    result.xyzw.x = a->xyzw.x * b->xyzw.x;
    result.xyzw.y = a->xyzw.y * b->xyzw.y;
    result.xyzw.z = a->xyzw.z * b->xyzw.z;
    result.xyzw.w = a->xyzw.w * b->xyzw.w;
    return result;
}

VECMATH_API fmat2 fmat2_mul(const fmat2* a, const fmat2* b)
{
    fmat2 result = { 0 };
    result.matrix.m00 = a->matrix.m00 * b->matrix.m00;
    result.matrix.m01 = a->matrix.m01 * b->matrix.m01;

    result.matrix.m10 = a->matrix.m10 * b->matrix.m10;
    result.matrix.m11 = a->matrix.m11 * b->matrix.m11;
    return result;
}

VECMATH_API fmat3 fmat3_mul(const fmat3* a, const fmat3* b)
{
    fmat3 result = { 0 };
    result.matrix.m00 = a->matrix.m00 * b->matrix.m00;
    result.matrix.m01 = a->matrix.m01 * b->matrix.m01;
    result.matrix.m02 = a->matrix.m02 * b->matrix.m02;

    result.matrix.m10 = a->matrix.m10 * b->matrix.m10;
    result.matrix.m11 = a->matrix.m11 * b->matrix.m11;
    result.matrix.m12 = a->matrix.m12 * b->matrix.m12;

    result.matrix.m20 = a->matrix.m20 * b->matrix.m20;
    result.matrix.m21 = a->matrix.m21 * b->matrix.m21;
    result.matrix.m22 = a->matrix.m22 * b->matrix.m22;
    return result;
}

VECMATH_API fmat4 fmat4_mul(const fmat4* a, const fmat4* b)
{
    fmat4 result = { 0 };
    result.matrix.m00 = a->matrix.m00 * b->matrix.m00;
    result.matrix.m01 = a->matrix.m01 * b->matrix.m01;
    result.matrix.m02 = a->matrix.m02 * b->matrix.m02;
    result.matrix.m03 = a->matrix.m03 * b->matrix.m03;

    result.matrix.m10 = a->matrix.m10 * b->matrix.m10;
    result.matrix.m11 = a->matrix.m11 * b->matrix.m11;
    result.matrix.m12 = a->matrix.m12 * b->matrix.m12;
    result.matrix.m13 = a->matrix.m13 * b->matrix.m13;

    result.matrix.m20 = a->matrix.m20 * b->matrix.m20;
    result.matrix.m21 = a->matrix.m21 * b->matrix.m21;
    result.matrix.m22 = a->matrix.m22 * b->matrix.m22;
    result.matrix.m23 = a->matrix.m23 * b->matrix.m23;

    result.matrix.m30 = a->matrix.m30 * b->matrix.m30;
    result.matrix.m31 = a->matrix.m31 * b->matrix.m31;
    result.matrix.m32 = a->matrix.m32 * b->matrix.m32;
    result.matrix.m33 = a->matrix.m33 * b->matrix.m33;
    return result;
}

VECMATH_API double2 double2_mul(const double2* a, const double2* b)
{
    double2 result = { 0 };
    result.xy.x = a->xy.x * b->xy.x;
    result.xy.y = a->xy.y * b->xy.y;
    return result;
}

double3 double3_mul(const double3* a, const double3* b)
{
    double3 result = { 0 };
    result.xyz.x = a->xyz.x * b->xyz.x;
    result.xyz.y = a->xyz.y * b->xyz.y;
    result.xyz.z = a->xyz.z * b->xyz.z;
    return result;
}

VECMATH_API double4 double4_mul(const double4* a, const double4* b)
{
    double4 result = { 0 };
    result.xyzw.x = a->xyzw.x * b->xyzw.x;
    result.xyzw.y = a->xyzw.y * b->xyzw.y;
    result.xyzw.z = a->xyzw.z * b->xyzw.z;
    result.xyzw.w = a->xyzw.w * b->xyzw.w;
    return result;
}

VECMATH_API dmat2 dmat2_mul(const dmat2* a, const dmat2* b)
{
    dmat2 result = { 0 };
    result.matrix.m00 = a->matrix.m00 * b->matrix.m00;
    result.matrix.m01 = a->matrix.m01 * b->matrix.m01;

    result.matrix.m10 = a->matrix.m10 * b->matrix.m10;
    result.matrix.m11 = a->matrix.m11 * b->matrix.m11;
    return result;
}

VECMATH_API dmat3 dmat3_mul(const dmat3* a, const dmat3* b)
{
    dmat3 result = { 0 };
    result.matrix.m00 = a->matrix.m00 * b->matrix.m00;
    result.matrix.m01 = a->matrix.m01 * b->matrix.m01;
    result.matrix.m02 = a->matrix.m02 * b->matrix.m02;

    result.matrix.m10 = a->matrix.m10 * b->matrix.m10;
    result.matrix.m11 = a->matrix.m11 * b->matrix.m11;
    result.matrix.m12 = a->matrix.m12 * b->matrix.m12;

    result.matrix.m20 = a->matrix.m20 * b->matrix.m20;
    result.matrix.m21 = a->matrix.m21 * b->matrix.m21;
    result.matrix.m22 = a->matrix.m22 * b->matrix.m22;
    return result;
}

VECMATH_API dmat4 dmat4_mul(const dmat4* a, const dmat4* b)
{
    dmat4 result = { 0 };
    result.matrix.m00 = a->matrix.m00 * b->matrix.m00;
    result.matrix.m01 = a->matrix.m01 * b->matrix.m01;
    result.matrix.m02 = a->matrix.m02 * b->matrix.m02;
    result.matrix.m03 = a->matrix.m03 * b->matrix.m03;

    result.matrix.m10 = a->matrix.m10 * b->matrix.m10;
    result.matrix.m11 = a->matrix.m11 * b->matrix.m11;
    result.matrix.m12 = a->matrix.m12 * b->matrix.m12;
    result.matrix.m13 = a->matrix.m13 * b->matrix.m13;

    result.matrix.m20 = a->matrix.m20 * b->matrix.m20;
    result.matrix.m21 = a->matrix.m21 * b->matrix.m21;
    result.matrix.m22 = a->matrix.m22 * b->matrix.m22;
    result.matrix.m23 = a->matrix.m23 * b->matrix.m23;

    result.matrix.m30 = a->matrix.m30 * b->matrix.m30;
    result.matrix.m31 = a->matrix.m31 * b->matrix.m31;
    result.matrix.m32 = a->matrix.m32 * b->matrix.m32;
    result.matrix.m33 = a->matrix.m33 * b->matrix.m33;
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// div
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API float2 float2_div(const float2* a, const float2* b)
{
    float2 result = { 0 };
    result.xy.x = a->xy.x / b->xy.x;
    result.xy.y = a->xy.y / b->xy.y;
    return result;
}

VECMATH_API float3 float3_div(const float3* a, const float3* b)
{
    float3 result = { 0 };
    result.xyz.x = a->xyz.x / b->xyz.x;
    result.xyz.y = a->xyz.y / b->xyz.y;
    result.xyz.z = a->xyz.z / b->xyz.z;
    return result;
}

VECMATH_API float4 float4_div(const float4* a, const float4* b)
{
    float4 result = { 0 };
    result.xyzw.x = a->xyzw.x * b->xyzw.x;
    result.xyzw.y = a->xyzw.y * b->xyzw.y;
    result.xyzw.z = a->xyzw.z * b->xyzw.z;
    result.xyzw.w = a->xyzw.w * b->xyzw.w;
    return result;
}

VECMATH_API fmat2 fmat2_div(const fmat2* a, const fmat2* b)
{
    fmat2 result = { 0 };
    result.matrix.m00 = a->matrix.m00 / b->matrix.m00;
    result.matrix.m01 = a->matrix.m01 / b->matrix.m01;

    result.matrix.m10 = a->matrix.m10 / b->matrix.m10;
    result.matrix.m11 = a->matrix.m11 / b->matrix.m11;
    return result;
}

VECMATH_API fmat3 fmat3_div(const fmat3* a, const fmat3* b)
{
    fmat3 result = { 0 };
    result.matrix.m00 = a->matrix.m00 / b->matrix.m00;
    result.matrix.m01 = a->matrix.m01 / b->matrix.m01;
    result.matrix.m02 = a->matrix.m02 / b->matrix.m02;

    result.matrix.m10 = a->matrix.m10 / b->matrix.m10;
    result.matrix.m11 = a->matrix.m11 / b->matrix.m11;
    result.matrix.m12 = a->matrix.m12 / b->matrix.m12;

    result.matrix.m20 = a->matrix.m20 / b->matrix.m20;
    result.matrix.m21 = a->matrix.m21 / b->matrix.m21;
    result.matrix.m22 = a->matrix.m22 / b->matrix.m22;
    return result;
}

VECMATH_API fmat4 fmat4_div(const fmat4* a, const fmat4* b)
{
    fmat4 result = { 0 };
    result.matrix.m00 = a->matrix.m00 / b->matrix.m00;
    result.matrix.m01 = a->matrix.m01 / b->matrix.m01;
    result.matrix.m02 = a->matrix.m02 / b->matrix.m02;
    result.matrix.m03 = a->matrix.m03 / b->matrix.m03;

    result.matrix.m10 = a->matrix.m10 / b->matrix.m10;
    result.matrix.m11 = a->matrix.m11 / b->matrix.m11;
    result.matrix.m12 = a->matrix.m12 / b->matrix.m12;
    result.matrix.m13 = a->matrix.m13 / b->matrix.m13;

    result.matrix.m20 = a->matrix.m20 / b->matrix.m20;
    result.matrix.m21 = a->matrix.m21 / b->matrix.m21;
    result.matrix.m22 = a->matrix.m22 / b->matrix.m22;
    result.matrix.m23 = a->matrix.m23 / b->matrix.m23;

    result.matrix.m30 = a->matrix.m30 / b->matrix.m30;
    result.matrix.m31 = a->matrix.m31 / b->matrix.m31;
    result.matrix.m32 = a->matrix.m32 / b->matrix.m32;
    result.matrix.m33 = a->matrix.m33 / b->matrix.m33;
    return result;
}

VECMATH_API double2 double2_div(const double2* a, const double2* b)
{
    double2 result = { 0 };
    result.xy.x = a->xy.x / b->xy.x;
    result.xy.y = a->xy.y / b->xy.y;
    return result;
}

VECMATH_API double3 double3_div(const double3* a, const double3* b)
{
    double3 result = { 0 };
    result.xyz.x = a->xyz.x / b->xyz.x;
    result.xyz.y = a->xyz.y / b->xyz.y;
    result.xyz.z = a->xyz.z / b->xyz.z;
    return result;
}

VECMATH_API double4 double4_div(const double4* a, const double4* b)
{
    double4 result = { 0 };
    result.xyzw.x = a->xyzw.x / b->xyzw.x;
    result.xyzw.y = a->xyzw.y / b->xyzw.y;
    result.xyzw.z = a->xyzw.z / b->xyzw.z;
    result.xyzw.w = a->xyzw.w / b->xyzw.w;
    return result;
}

VECMATH_API dmat2 dmat2_div(const dmat2* a, const dmat2* b)
{
    dmat2 result = { 0 };
    result.matrix.m00 = a->matrix.m00 / b->matrix.m00;
    result.matrix.m01 = a->matrix.m01 / b->matrix.m01;

    result.matrix.m10 = a->matrix.m10 / b->matrix.m10;
    result.matrix.m11 = a->matrix.m11 / b->matrix.m11;
    return result;
}

VECMATH_API dmat3 dmat3_div(const dmat3* a, const dmat3* b)
{
    dmat3 result = { 0 };
    result.matrix.m00 = a->matrix.m00 / b->matrix.m00;
    result.matrix.m01 = a->matrix.m01 / b->matrix.m01;
    result.matrix.m02 = a->matrix.m02 / b->matrix.m02;

    result.matrix.m10 = a->matrix.m10 / b->matrix.m10;
    result.matrix.m11 = a->matrix.m11 / b->matrix.m11;
    result.matrix.m12 = a->matrix.m12 / b->matrix.m12;

    result.matrix.m20 = a->matrix.m20 / b->matrix.m20;
    result.matrix.m21 = a->matrix.m21 / b->matrix.m21;
    result.matrix.m22 = a->matrix.m22 / b->matrix.m22;
    return result;
}

VECMATH_API dmat4 dmat4_div(const dmat4* a, const dmat4* b)
{
    dmat4 result = { 0 };
    result.matrix.m00 = a->matrix.m00 / b->matrix.m00;
    result.matrix.m01 = a->matrix.m01 / b->matrix.m01;
    result.matrix.m02 = a->matrix.m02 / b->matrix.m02;
    result.matrix.m03 = a->matrix.m03 / b->matrix.m03;

    result.matrix.m10 = a->matrix.m10 / b->matrix.m10;
    result.matrix.m11 = a->matrix.m11 / b->matrix.m11;
    result.matrix.m12 = a->matrix.m12 / b->matrix.m12;
    result.matrix.m13 = a->matrix.m13 / b->matrix.m13;

    result.matrix.m20 = a->matrix.m20 / b->matrix.m20;
    result.matrix.m21 = a->matrix.m21 / b->matrix.m21;
    result.matrix.m22 = a->matrix.m22 / b->matrix.m22;
    result.matrix.m23 = a->matrix.m23 / b->matrix.m23;

    result.matrix.m30 = a->matrix.m30 / b->matrix.m30;
    result.matrix.m31 = a->matrix.m31 / b->matrix.m31;
    result.matrix.m32 = a->matrix.m32 / b->matrix.m32;
    result.matrix.m33 = a->matrix.m33 / b->matrix.m33;
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// equals
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API vecbool float2_equals(const float2 *a, const float2 *b)
{
    return memcmp(a->data, b->data, sizeof(a->data)) == 0;
}

VECMATH_API vecbool float3_equals(const float3 *a, const float3 *b)
{
    return memcmp(a->data, b->data, sizeof(a->data)) == 0;
}

VECMATH_API vecbool float4_equals(const float4 *a, const float4 *b)
{
    return memcmp(a->data, b->data, sizeof(a->data)) == 0;
}

VECMATH_API vecbool fmat2_equals(const fmat2 *a, const fmat2 *b)
{
    return memcmp(a->data, b->data, sizeof(a->data)) == 0;
}

VECMATH_API vecbool fmat3_equals(const fmat3 *a, const fmat3 *b)
{
    return memcmp(a->data, b->data, sizeof(a->data)) == 0;
}

VECMATH_API vecbool fmat4_equals(const fmat4 *a, const fmat4 *b)
{
    return memcmp(a->data, b->data, sizeof(a->data)) == 0;
}

VECMATH_API vecbool double2_equals(const double2 *a, const double2 *b)
{
    return memcmp(a->data, b->data, sizeof(a->data)) == 0;
}

VECMATH_API vecbool double3_equals(const double3 *a, const double3 *b)
{
    return memcmp(a->data, b->data, sizeof(a->data)) == 0;
}

VECMATH_API vecbool double4_equals(const double4 *a, const double4 *b)
{
    return memcmp(a->data, b->data, sizeof(a->data)) == 0;
}

VECMATH_API vecbool dmat2_equals(const dmat2 *a, const dmat2 *b)
{
    return memcmp(a->data, b->data, sizeof(a->data)) == 0;
}

VECMATH_API vecbool dmat3_equals(const dmat3 *a, const dmat3 *b)
{
    return memcmp(a->data, b->data, sizeof(a->data)) == 0;
}

VECMATH_API vecbool dmat4_equals(const dmat4 *a, const dmat4 *b)
{
    return memcmp(a->data, b->data, sizeof(a->data)) == 0;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// aprox_equals
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API vecbool float2_aprox_equals(const float2 *a, const float2 *b)
{
    return  (fabsf(a->xy.x - b->xy.x) <= VECMATH_EPSILON_FZERO) &&
            (fabsf(a->xy.y - b->xy.y) <= VECMATH_EPSILON_FZERO);
}

VECMATH_API vecbool float3_aprox_equals(const float3 *a, const float3 *b)
{
    return  (fabsf(a->xyz.x - b->xyz.x) <= VECMATH_EPSILON_FZERO) &&
            (fabsf(a->xyz.y - b->xyz.y) <= VECMATH_EPSILON_FZERO) && 
            (fabsf(a->xyz.z - b->xyz.z) <= VECMATH_EPSILON_FZERO);
}

VECMATH_API vecbool float4_aprox_equals(const float4 *a, const float4 *b)
{
    return  (fabsf(a->xyzw.x - b->xyzw.x) <= VECMATH_EPSILON_FZERO) &&
            (fabsf(a->xyzw.y - b->xyzw.y) <= VECMATH_EPSILON_FZERO) && 
            (fabsf(a->xyzw.z - b->xyzw.z) <= VECMATH_EPSILON_FZERO) &&
            (fabsf(a->xyzw.w - b->xyzw.w) <= VECMATH_EPSILON_FZERO);
}

VECMATH_API vecbool fmat2_aprox_equals(const fmat2 *a, const fmat2 *b)
{
    return  (fabsf(a->matrix.m00 - b->matrix.m00) <= VECMATH_EPSILON_FZERO) &&
            (fabsf(a->matrix.m01 - b->matrix.m01) <= VECMATH_EPSILON_FZERO) && 

            (fabsf(a->matrix.m10 - b->matrix.m10) <= VECMATH_EPSILON_FZERO) &&
            (fabsf(a->matrix.m11 - b->matrix.m11) <= VECMATH_EPSILON_FZERO);
}

VECMATH_API vecbool fmat3_aprox_equals(const fmat3 *a, const fmat3 *b)
{
    return  (fabsf(a->matrix.m00 - b->matrix.m00) <= VECMATH_EPSILON_FZERO) &&
            (fabsf(a->matrix.m01 - b->matrix.m01) <= VECMATH_EPSILON_FZERO) && 
            (fabsf(a->matrix.m02 - b->matrix.m02) <= VECMATH_EPSILON_FZERO) && 

            (fabsf(a->matrix.m10 - b->matrix.m10) <= VECMATH_EPSILON_FZERO) &&
            (fabsf(a->matrix.m11 - b->matrix.m11) <= VECMATH_EPSILON_FZERO) &&
            (fabsf(a->matrix.m12 - b->matrix.m12) <= VECMATH_EPSILON_FZERO) && 

            (fabsf(a->matrix.m20 - b->matrix.m20) <= VECMATH_EPSILON_FZERO) &&
            (fabsf(a->matrix.m21 - b->matrix.m21) <= VECMATH_EPSILON_FZERO) &&
            (fabsf(a->matrix.m22 - b->matrix.m22) <= VECMATH_EPSILON_FZERO);
}

VECMATH_API vecbool fmat4_aprox_equals(const fmat4 *a, const fmat4 *b)
{
    return  (fabsf(a->matrix.m00 - b->matrix.m00) <= VECMATH_EPSILON_FZERO) &&
            (fabsf(a->matrix.m01 - b->matrix.m01) <= VECMATH_EPSILON_FZERO) && 
            (fabsf(a->matrix.m02 - b->matrix.m02) <= VECMATH_EPSILON_FZERO) && 
            (fabsf(a->matrix.m03 - b->matrix.m03) <= VECMATH_EPSILON_FZERO) && 

            (fabsf(a->matrix.m10 - b->matrix.m10) <= VECMATH_EPSILON_FZERO) &&
            (fabsf(a->matrix.m11 - b->matrix.m11) <= VECMATH_EPSILON_FZERO) &&
            (fabsf(a->matrix.m12 - b->matrix.m12) <= VECMATH_EPSILON_FZERO) && 
            (fabsf(a->matrix.m13 - b->matrix.m13) <= VECMATH_EPSILON_FZERO) && 

            (fabsf(a->matrix.m20 - b->matrix.m20) <= VECMATH_EPSILON_FZERO) &&
            (fabsf(a->matrix.m21 - b->matrix.m21) <= VECMATH_EPSILON_FZERO) &&
            (fabsf(a->matrix.m22 - b->matrix.m22) <= VECMATH_EPSILON_FZERO) &&
            (fabsf(a->matrix.m23 - b->matrix.m23) <= VECMATH_EPSILON_FZERO) && 

            (fabsf(a->matrix.m30 - b->matrix.m30) <= VECMATH_EPSILON_FZERO) &&
            (fabsf(a->matrix.m31 - b->matrix.m31) <= VECMATH_EPSILON_FZERO) &&
            (fabsf(a->matrix.m32 - b->matrix.m32) <= VECMATH_EPSILON_FZERO) &&
            (fabsf(a->matrix.m33 - b->matrix.m33) <= VECMATH_EPSILON_FZERO);
}

VECMATH_API vecbool double2_aprox_equals(const double2 *a, const double2 *b)
{
    return  (fabs(a->xy.x - b->xy.x) <= VECMATH_EPSILON_DZERO) &&
            (fabs(a->xy.y - b->xy.y) <= VECMATH_EPSILON_DZERO);
}

VECMATH_API vecbool double3_aprox_equals(const double3 *a, const double3 *b)
{
    return  (fabs(a->xyz.x - b->xyz.x) <= VECMATH_EPSILON_DZERO) &&
            (fabs(a->xyz.y - b->xyz.y) <= VECMATH_EPSILON_DZERO) && 
            (fabs(a->xyz.z - b->xyz.z) <= VECMATH_EPSILON_DZERO);
}

VECMATH_API vecbool double4_aprox_equals(const double4 *a, const double4 *b)
{
    return  (fabs(a->xyzw.x - b->xyzw.x) <= VECMATH_EPSILON_DZERO) &&
            (fabs(a->xyzw.y - b->xyzw.y) <= VECMATH_EPSILON_DZERO) && 
            (fabs(a->xyzw.z - b->xyzw.z) <= VECMATH_EPSILON_DZERO) &&
            (fabs(a->xyzw.w - b->xyzw.w) <= VECMATH_EPSILON_DZERO);
}

VECMATH_API vecbool dmat2_aprox_equals(const dmat2 *a, const dmat2 *b)
{
    return  (fabs(a->matrix.m00 - b->matrix.m00) <= VECMATH_EPSILON_DZERO) &&
            (fabs(a->matrix.m01 - b->matrix.m01) <= VECMATH_EPSILON_DZERO) && 

            (fabs(a->matrix.m10 - b->matrix.m10) <= VECMATH_EPSILON_DZERO) &&
            (fabs(a->matrix.m11 - b->matrix.m11) <= VECMATH_EPSILON_DZERO);
}

VECMATH_API vecbool dmat3_aprox_equals(const dmat3 *a, const dmat3 *b)
{
    return  (fabs(a->matrix.m00 - b->matrix.m00) <= VECMATH_EPSILON_DZERO) &&
            (fabs(a->matrix.m01 - b->matrix.m01) <= VECMATH_EPSILON_DZERO) && 
            (fabs(a->matrix.m02 - b->matrix.m02) <= VECMATH_EPSILON_DZERO) && 

            (fabs(a->matrix.m10 - b->matrix.m10) <= VECMATH_EPSILON_DZERO) &&
            (fabs(a->matrix.m11 - b->matrix.m11) <= VECMATH_EPSILON_DZERO) &&
            (fabs(a->matrix.m12 - b->matrix.m12) <= VECMATH_EPSILON_DZERO) && 

            (fabs(a->matrix.m20 - b->matrix.m20) <= VECMATH_EPSILON_DZERO) &&
            (fabs(a->matrix.m21 - b->matrix.m21) <= VECMATH_EPSILON_DZERO) &&
            (fabs(a->matrix.m22 - b->matrix.m22) <= VECMATH_EPSILON_DZERO);
}

VECMATH_API vecbool dmat4_aprox_equals(const dmat4 *a, const dmat4 *b)
{
    return  (fabs(a->matrix.m00 - b->matrix.m00) <= VECMATH_EPSILON_DZERO) &&
            (fabs(a->matrix.m01 - b->matrix.m01) <= VECMATH_EPSILON_DZERO) && 
            (fabs(a->matrix.m02 - b->matrix.m02) <= VECMATH_EPSILON_DZERO) && 
            (fabs(a->matrix.m03 - b->matrix.m03) <= VECMATH_EPSILON_DZERO) && 

            (fabs(a->matrix.m10 - b->matrix.m10) <= VECMATH_EPSILON_DZERO) &&
            (fabs(a->matrix.m11 - b->matrix.m11) <= VECMATH_EPSILON_DZERO) &&
            (fabs(a->matrix.m12 - b->matrix.m12) <= VECMATH_EPSILON_DZERO) && 
            (fabs(a->matrix.m13 - b->matrix.m13) <= VECMATH_EPSILON_DZERO) && 

            (fabs(a->matrix.m20 - b->matrix.m20) <= VECMATH_EPSILON_DZERO) &&
            (fabs(a->matrix.m21 - b->matrix.m21) <= VECMATH_EPSILON_DZERO) &&
            (fabs(a->matrix.m22 - b->matrix.m22) <= VECMATH_EPSILON_DZERO) &&
            (fabs(a->matrix.m23 - b->matrix.m23) <= VECMATH_EPSILON_DZERO) && 

            (fabs(a->matrix.m30 - b->matrix.m30) <= VECMATH_EPSILON_DZERO) &&
            (fabs(a->matrix.m31 - b->matrix.m31) <= VECMATH_EPSILON_DZERO) &&
            (fabs(a->matrix.m32 - b->matrix.m32) <= VECMATH_EPSILON_DZERO) &&
            (fabs(a->matrix.m33 - b->matrix.m33) <= VECMATH_EPSILON_DZERO);
}

