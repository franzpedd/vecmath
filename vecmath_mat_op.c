#include "vecmath_mat_op.h"
#include "vecmath_basic_op.h"
#include "vecmath_vec_op.h"
#include <math.h>

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// identity 
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API fmat2 fmat2_identity()
{
    fmat2 mat;
    mat.matrix.m00 = 1.0f; mat.matrix.m01 = 0.0f;
    mat.matrix.m10 = 0.0f; mat.matrix.m11 = 1.0f;
    return mat;
}

VECMATH_API fmat3 fmat3_identity()
{
    fmat3 mat;
    mat.matrix.m00 = 1.0f; mat.matrix.m01 = 0.0f; mat.matrix.m02 = 0.0f;
    mat.matrix.m10 = 0.0f; mat.matrix.m11 = 1.0f; mat.matrix.m12 = 0.0f;
    mat.matrix.m20 = 0.0f; mat.matrix.m21 = 0.0f; mat.matrix.m22 = 1.0f;
    return mat;
}

VECMATH_API fmat4 fmat4_identity()
{
    fmat4 mat;
    mat.matrix.m00 = 1.0f; mat.matrix.m01 = 0.0f; mat.matrix.m02 = 0.0f; mat.matrix.m03 = 0.0f;
    mat.matrix.m10 = 0.0f; mat.matrix.m11 = 1.0f; mat.matrix.m12 = 0.0f; mat.matrix.m13 = 0.0f;
    mat.matrix.m20 = 0.0f; mat.matrix.m21 = 0.0f; mat.matrix.m22 = 1.0f; mat.matrix.m23 = 0.0f;
    mat.matrix.m30 = 0.0f; mat.matrix.m31 = 0.0f; mat.matrix.m32 = 0.0f; mat.matrix.m33 = 1.0f;
    return mat;
}

VECMATH_API dmat2 dmat2_identity()
{
    dmat2 mat;
    mat.matrix.m00 = 1.0f; mat.matrix.m01 = 0.0f;
    mat.matrix.m10 = 0.0f; mat.matrix.m11 = 1.0f;
    return mat;
}

VECMATH_API dmat3 dmat3_identity()
{
    dmat3 mat;
    mat.matrix.m00 = 1.0f; mat.matrix.m01 = 0.0f; mat.matrix.m02 = 0.0f;
    mat.matrix.m10 = 0.0f; mat.matrix.m11 = 1.0f; mat.matrix.m12 = 0.0f;
    mat.matrix.m20 = 0.0f; mat.matrix.m21 = 0.0f; mat.matrix.m22 = 1.0f;
    return mat;
}

VECMATH_API dmat4 dmat4_identity()
{
    dmat4 mat;
    mat.matrix.m00 = 1.0f; mat.matrix.m01 = 0.0f; mat.matrix.m02 = 0.0f; mat.matrix.m03 = 0.0f;
    mat.matrix.m10 = 0.0f; mat.matrix.m11 = 1.0f; mat.matrix.m12 = 0.0f; mat.matrix.m13 = 0.0f;
    mat.matrix.m20 = 0.0f; mat.matrix.m21 = 0.0f; mat.matrix.m22 = 1.0f; mat.matrix.m23 = 0.0f;
    mat.matrix.m30 = 0.0f; mat.matrix.m31 = 0.0f; mat.matrix.m32 = 0.0f; mat.matrix.m33 = 1.0f;
    return mat;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// mul_vec 
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API fmat2 fmat2_mul_float2(const fmat2* m, const float2* v)
{
    fmat2 result = { 0 };
    result.matrix.m00 = m->matrix.m00 * v->xy.x;
    result.matrix.m10 = m->matrix.m10 * v->xy.x;
    
    result.matrix.m01 = m->matrix.m01 * v->xy.y;
    result.matrix.m11 = m->matrix.m11 * v->xy.y;
    return result;
}

VECMATH_API fmat3 fmat3_mul_float3(const fmat3* m, const float3* v)
{
    fmat3 result = { 0 };
    result.matrix.m00 = m->matrix.m00 * v->xyz.x;
    result.matrix.m10 = m->matrix.m10 * v->xyz.x;
    result.matrix.m20 = m->matrix.m20 * v->xyz.x;
    
    result.matrix.m01 = m->matrix.m01 * v->xyz.y;
    result.matrix.m11 = m->matrix.m11 * v->xyz.y;
    result.matrix.m21 = m->matrix.m21 * v->xyz.y;
    
    result.matrix.m02 = m->matrix.m02 * v->xyz.z;
    result.matrix.m12 = m->matrix.m12 * v->xyz.z;
    result.matrix.m22 = m->matrix.m22 * v->xyz.z;
    return result;
}

VECMATH_API fmat4 fmat4_mul_float4(const fmat4* m, const float4* v)
{
    fmat4 result = { 0 };
    result.matrix.m00 = m->matrix.m00 * v->xyzw.x;
    result.matrix.m10 = m->matrix.m10 * v->xyzw.x;
    result.matrix.m20 = m->matrix.m20 * v->xyzw.x;
    result.matrix.m30 = m->matrix.m30 * v->xyzw.x;
    
    result.matrix.m01 = m->matrix.m01 * v->xyzw.y;
    result.matrix.m11 = m->matrix.m11 * v->xyzw.y;
    result.matrix.m21 = m->matrix.m21 * v->xyzw.y;
    result.matrix.m31 = m->matrix.m31 * v->xyzw.y;
    
    result.matrix.m02 = m->matrix.m02 * v->xyzw.z;
    result.matrix.m12 = m->matrix.m12 * v->xyzw.z;
    result.matrix.m22 = m->matrix.m22 * v->xyzw.z;
    result.matrix.m32 = m->matrix.m32 * v->xyzw.z;
    
    result.matrix.m03 = m->matrix.m03 * v->xyzw.w;
    result.matrix.m13 = m->matrix.m13 * v->xyzw.w;
    result.matrix.m23 = m->matrix.m23 * v->xyzw.w;
    result.matrix.m33 = m->matrix.m33 * v->xyzw.w;
    return result;
}

VECMATH_API dmat2 dmat2_mul_double2(const dmat2* m, const double2* v)
{
    dmat2 result = { 0 };
    result.matrix.m00 = m->matrix.m00 * v->xy.x;
    result.matrix.m10 = m->matrix.m10 * v->xy.x;
    
    result.matrix.m01 = m->matrix.m01 * v->xy.y;
    result.matrix.m11 = m->matrix.m11 * v->xy.y;
    return result;
}

VECMATH_API dmat3 dmat3_mul_double3(const dmat3* m, const double3* v)
{
    dmat3 result = { 0 };
    result.matrix.m00 = m->matrix.m00 * v->xyz.x;
    result.matrix.m10 = m->matrix.m10 * v->xyz.x;
    result.matrix.m20 = m->matrix.m20 * v->xyz.x;
    
    result.matrix.m01 = m->matrix.m01 * v->xyz.y;
    result.matrix.m11 = m->matrix.m11 * v->xyz.y;
    result.matrix.m21 = m->matrix.m21 * v->xyz.y;
    
    result.matrix.m02 = m->matrix.m02 * v->xyz.z;
    result.matrix.m12 = m->matrix.m12 * v->xyz.z;
    result.matrix.m22 = m->matrix.m22 * v->xyz.z;
    return result;
}

VECMATH_API dmat4 dmat4_mul_double4(const dmat4* m, const double4* v)
{
    dmat4 result = { 0 };
    result.matrix.m00 = m->matrix.m00 * v->xyzw.x;
    result.matrix.m10 = m->matrix.m10 * v->xyzw.x;
    result.matrix.m20 = m->matrix.m20 * v->xyzw.x;
    result.matrix.m30 = m->matrix.m30 * v->xyzw.x;
    
    result.matrix.m01 = m->matrix.m01 * v->xyzw.y;
    result.matrix.m11 = m->matrix.m11 * v->xyzw.y;
    result.matrix.m21 = m->matrix.m21 * v->xyzw.y;
    result.matrix.m31 = m->matrix.m31 * v->xyzw.y;
    
    result.matrix.m02 = m->matrix.m02 * v->xyzw.z;
    result.matrix.m12 = m->matrix.m12 * v->xyzw.z;
    result.matrix.m22 = m->matrix.m22 * v->xyzw.z;
    result.matrix.m32 = m->matrix.m32 * v->xyzw.z;
    
    result.matrix.m03 = m->matrix.m03 * v->xyzw.w;
    result.matrix.m13 = m->matrix.m13 * v->xyzw.w;
    result.matrix.m23 = m->matrix.m23 * v->xyzw.w;
    result.matrix.m33 = m->matrix.m33 * v->xyzw.w;
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// transpose 
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API fmat2 fmat2_transpose(const fmat2* m)
{
    fmat2 result = { 0 };
    result.matrix.m00 = m->matrix.m00; 
    result.matrix.m01 = m->matrix.m10;

    result.matrix.m10 = m->matrix.m01;
    result.matrix.m11 = m->matrix.m11;
    return result;
}

VECMATH_API fmat3 fmat3_transpose(const fmat3* m)
{
    fmat3 result = { 0 };
    result.matrix.m00 = m->matrix.m00;
    result.matrix.m01 = m->matrix.m10;
    result.matrix.m02 = m->matrix.m20;
    
    result.matrix.m10 = m->matrix.m01;
    result.matrix.m11 = m->matrix.m11;
    result.matrix.m12 = m->matrix.m21;
    
    result.matrix.m20 = m->matrix.m02;
    result.matrix.m21 = m->matrix.m12;
    result.matrix.m22 = m->matrix.m22;
    return result;
}

VECMATH_API fmat4 fmat4_transpose(const fmat4* m)
{
    fmat4 result = { 0 };
    result.matrix.m00 = m->matrix.m00;
    result.matrix.m01 = m->matrix.m10;
    result.matrix.m02 = m->matrix.m20;
    result.matrix.m03 = m->matrix.m30;
    
    result.matrix.m10 = m->matrix.m01;
    result.matrix.m11 = m->matrix.m11;
    result.matrix.m12 = m->matrix.m21;
    result.matrix.m13 = m->matrix.m31;
    
    result.matrix.m20 = m->matrix.m02;
    result.matrix.m21 = m->matrix.m12;
    result.matrix.m22 = m->matrix.m22;
    result.matrix.m23 = m->matrix.m32;
    
    result.matrix.m30 = m->matrix.m03;
    result.matrix.m31 = m->matrix.m13;
    result.matrix.m32 = m->matrix.m23;
    result.matrix.m33 = m->matrix.m33;
    return result;
}

VECMATH_API dmat2 dmat2_transpose(const dmat2* m)
{
    dmat2 result = { 0 };
    result.matrix.m00 = m->matrix.m00; 
    result.matrix.m01 = m->matrix.m10;

    result.matrix.m10 = m->matrix.m01;
    result.matrix.m11 = m->matrix.m11;
    return result;
}

VECMATH_API dmat3 dmat3_transpose(const dmat3* m)
{
    dmat3 result = { 0 };
    result.matrix.m00 = m->matrix.m00;
    result.matrix.m01 = m->matrix.m10;
    result.matrix.m02 = m->matrix.m20;
    
    result.matrix.m10 = m->matrix.m01;
    result.matrix.m11 = m->matrix.m11;
    result.matrix.m12 = m->matrix.m21;
    
    result.matrix.m20 = m->matrix.m02;
    result.matrix.m21 = m->matrix.m12;
    result.matrix.m22 = m->matrix.m22;
    return result;
}

VECMATH_API dmat4 dmat4_transpose(const dmat4* m)
{
    dmat4 result = { 0 };
    result.matrix.m00 = m->matrix.m00;
    result.matrix.m01 = m->matrix.m10;
    result.matrix.m02 = m->matrix.m20;
    result.matrix.m03 = m->matrix.m30;
    
    result.matrix.m10 = m->matrix.m01;
    result.matrix.m11 = m->matrix.m11;
    result.matrix.m12 = m->matrix.m21;
    result.matrix.m13 = m->matrix.m31;
    
    result.matrix.m20 = m->matrix.m02;
    result.matrix.m21 = m->matrix.m12;
    result.matrix.m22 = m->matrix.m22;
    result.matrix.m23 = m->matrix.m32;
    
    result.matrix.m30 = m->matrix.m03;
    result.matrix.m31 = m->matrix.m13;
    result.matrix.m32 = m->matrix.m23;
    result.matrix.m33 = m->matrix.m33;
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// determinant 
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API float fmat2_determinant(const fmat2* m)
{
    return (m->matrix.m00 * m->matrix.m11) - (m->matrix.m01 * m->matrix.m10);
}

VECMATH_API float fmat3_determinant(const fmat3* m)
{
    return m->matrix.m00 * (m->matrix.m11 * m->matrix.m22 - m->matrix.m12 * m->matrix.m21) 
        - m->matrix.m01 * (m->matrix.m10 * m->matrix.m22 - m->matrix.m12 * m->matrix.m20) 
        + m->matrix.m02 * (m->matrix.m10 * m->matrix.m21 - m->matrix.m11 * m->matrix.m20);
}

VECMATH_API float fmat4_determinant(const fmat4* m)
{
    float det = 0;
    float sub[3][3];
    
    for (int x = 0; x < 4; x++) {
        int subi = 0;
        for (int i = 1; i < 4; i++) {
            int subj = 0;
            for (int j = 0; j < 4; j++) {
                if (j == x) continue;
                sub[subi][subj] = m->data[i][j];
                subj++;
            }
            subi++;
        }
        fmat3 submat = {{{sub[0][0], sub[0][1], sub[0][2]},
            {sub[1][0], sub[1][1], sub[1][2]},
            {sub[2][0], sub[2][1], sub[2][2]}}};
        det += (x % 2 == 0 ? 1 : -1) * m->data[0][x] * fmat3_determinant(&submat);
    }
    return det;
}

VECMATH_API double dmat2_determinant(const dmat2* m)
{
    return (m->matrix.m00 * m->matrix.m11) - (m->matrix.m01 * m->matrix.m10);
}

VECMATH_API double dmat3_determinant(const dmat3* m)
{
    return m->matrix.m00 * (m->matrix.m11 * m->matrix.m22 - m->matrix.m12 * m->matrix.m21) 
        - m->matrix.m01 * (m->matrix.m10 * m->matrix.m22 - m->matrix.m12 * m->matrix.m20) 
        + m->matrix.m02 * (m->matrix.m10 * m->matrix.m21 - m->matrix.m11 * m->matrix.m20);
}

VECMATH_API double dmat4_determinant(const dmat4* m)
{
    double det = 0;
    double sub[3][3];
    
    for (int x = 0; x < 4; x++) {
        int subi = 0;
        for (int i = 1; i < 4; i++) {
            int subj = 0;
            for (int j = 0; j < 4; j++) {
                if (j == x) continue;
                sub[subi][subj] = m->data[i][j];
                subj++;
            }
            subi++;
        }
        dmat3 submat = {{{sub[0][0], sub[0][1], sub[0][2]},
            {sub[1][0], sub[1][1], sub[1][2]},
            {sub[2][0], sub[2][1], sub[2][2]}}};
        det += (x % 2 == 0 ? 1 : -1) * m->data[0][x] * dmat3_determinant(&submat);
    }
    return det;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// inverse 
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API fmat2 fmat2_inverse(const fmat2* m)
{
    fmat2 result = { 0 };
    float det = fmat2_determinant(m);

    if (fabsf(det) < VECMATH_EPSILON_FZERO) {
        return fmat2_identity();
    }
    
    float inv_det = 1.0f / det;
    result.matrix.m00 =  m->matrix.m11 * inv_det;
    result.matrix.m01 = -m->matrix.m01 * inv_det;
    result.matrix.m10 = -m->matrix.m10 * inv_det;
    result.matrix.m11 =  m->matrix.m00 * inv_det;
    return result;
}

VECMATH_API fmat3 fmat3_inverse(const fmat3* m)
{
    fmat3 result = { 0 };
    float det = fmat3_determinant(m);
    if (fabsf(det) < VECMATH_EPSILON_FZERO) {
        return fmat3_identity();
    }
    
    float inv_det = 1.0f / det;
    result.matrix.m00 = (m->matrix.m11 * m->matrix.m22 - m->matrix.m12 * m->matrix.m21) * inv_det;
    result.matrix.m01 = (m->matrix.m02 * m->matrix.m21 - m->matrix.m01 * m->matrix.m22) * inv_det;
    result.matrix.m02 = (m->matrix.m01 * m->matrix.m12 - m->matrix.m02 * m->matrix.m11) * inv_det;
    
    result.matrix.m10 = (m->matrix.m12 * m->matrix.m20 - m->matrix.m10 * m->matrix.m22) * inv_det;
    result.matrix.m11 = (m->matrix.m00 * m->matrix.m22 - m->matrix.m02 * m->matrix.m20) * inv_det;
    result.matrix.m12 = (m->matrix.m02 * m->matrix.m10 - m->matrix.m00 * m->matrix.m12) * inv_det;
    
    result.matrix.m20 = (m->matrix.m10 * m->matrix.m21 - m->matrix.m11 * m->matrix.m20) * inv_det;
    result.matrix.m21 = (m->matrix.m01 * m->matrix.m20 - m->matrix.m00 * m->matrix.m21) * inv_det;
    result.matrix.m22 = (m->matrix.m00 * m->matrix.m11 - m->matrix.m01 * m->matrix.m10) * inv_det;
    return result;
}

VECMATH_API fmat4 fmat4_inverse(const fmat4* m)
{
    // Using your matrix struct for clarity
    const float* mm = &m->matrix.m00;
    
    float inv[16];
    
    inv[0] = mm[5]  * mm[10] * mm[15] - 
             mm[5]  * mm[11] * mm[14] - 
             mm[9]  * mm[6]  * mm[15] + 
             mm[9]  * mm[7]  * mm[14] +
             mm[13] * mm[6]  * mm[11] - 
             mm[13] * mm[7]  * mm[10];

    inv[4] = -mm[4]  * mm[10] * mm[15] + 
              mm[4]  * mm[11] * mm[14] + 
              mm[8]  * mm[6]  * mm[15] - 
              mm[8]  * mm[7]  * mm[14] - 
              mm[12] * mm[6]  * mm[11] + 
              mm[12] * mm[7]  * mm[10];

    inv[8] = mm[4]  * mm[9] * mm[15] - 
             mm[4]  * mm[11] * mm[13] - 
             mm[8]  * mm[5] * mm[15] + 
             mm[8]  * mm[7] * mm[13] + 
             mm[12] * mm[5] * mm[11] - 
             mm[12] * mm[7] * mm[9];

    inv[12] = -mm[4]  * mm[9] * mm[14] + 
               mm[4]  * mm[10] * mm[13] +
               mm[8]  * mm[5] * mm[14] - 
               mm[8]  * mm[6] * mm[13] - 
               mm[12] * mm[5] * mm[10] + 
               mm[12] * mm[6] * mm[9];

    inv[1] = -mm[1]  * mm[10] * mm[15] + 
              mm[1]  * mm[11] * mm[14] + 
              mm[9]  * mm[2] * mm[15] - 
              mm[9]  * mm[3] * mm[14] - 
              mm[13] * mm[2] * mm[11] + 
              mm[13] * mm[3] * mm[10];

    inv[5] = mm[0]  * mm[10] * mm[15] - 
             mm[0]  * mm[11] * mm[14] - 
             mm[8]  * mm[2] * mm[15] + 
             mm[8]  * mm[3] * mm[14] + 
             mm[12] * mm[2] * mm[11] - 
             mm[12] * mm[3] * mm[10];

    inv[9] = -mm[0]  * mm[9] * mm[15] + 
              mm[0]  * mm[11] * mm[13] + 
              mm[8]  * mm[1] * mm[15] - 
              mm[8]  * mm[3] * mm[13] - 
              mm[12] * mm[1] * mm[11] + 
              mm[12] * mm[3] * mm[9];

    inv[13] = mm[0]  * mm[9] * mm[14] - 
              mm[0]  * mm[10] * mm[13] - 
              mm[8]  * mm[1] * mm[14] + 
              mm[8]  * mm[2] * mm[13] + 
              mm[12] * mm[1] * mm[10] - 
              mm[12] * mm[2] * mm[9];

    inv[2] = mm[1]  * mm[6] * mm[15] - 
             mm[1]  * mm[7] * mm[14] - 
             mm[5]  * mm[2] * mm[15] + 
             mm[5]  * mm[3] * mm[14] + 
             mm[13] * mm[2] * mm[7] - 
             mm[13] * mm[3] * mm[6];

    inv[6] = -mm[0]  * mm[6] * mm[15] + 
              mm[0]  * mm[7] * mm[14] + 
              mm[4]  * mm[2] * mm[15] - 
              mm[4]  * mm[3] * mm[14] - 
              mm[12] * mm[2] * mm[7] + 
              mm[12] * mm[3] * mm[6];

    inv[10] = mm[0]  * mm[5] * mm[15] - 
              mm[0]  * mm[7] * mm[13] - 
              mm[4]  * mm[1] * mm[15] + 
              mm[4]  * mm[3] * mm[13] + 
              mm[12] * mm[1] * mm[7] - 
              mm[12] * mm[3] * mm[5];

    inv[14] = -mm[0]  * mm[5] * mm[14] + 
               mm[0]  * mm[6] * mm[13] + 
               mm[4]  * mm[1] * mm[14] - 
               mm[4]  * mm[2] * mm[13] - 
               mm[12] * mm[1] * mm[6] + 
               mm[12] * mm[2] * mm[5];

    inv[3] = -mm[1] * mm[6] * mm[11] + 
              mm[1] * mm[7] * mm[10] + 
              mm[5] * mm[2] * mm[11] - 
              mm[5] * mm[3] * mm[10] - 
              mm[9] * mm[2] * mm[7] + 
              mm[9] * mm[3] * mm[6];

    inv[7] = mm[0] * mm[6] * mm[11] - 
             mm[0] * mm[7] * mm[10] - 
             mm[4] * mm[2] * mm[11] + 
             mm[4] * mm[3] * mm[10] + 
             mm[8] * mm[2] * mm[7] - 
             mm[8] * mm[3] * mm[6];

    inv[11] = -mm[0] * mm[5] * mm[11] + 
               mm[0] * mm[7] * mm[9] + 
               mm[4] * mm[1] * mm[11] - 
               mm[4] * mm[3] * mm[9] - 
               mm[8] * mm[1] * mm[7] + 
               mm[8] * mm[3] * mm[5];

    inv[15] = mm[0] * mm[5] * mm[10] - 
              mm[0] * mm[6] * mm[9] - 
              mm[4] * mm[1] * mm[10] + 
              mm[4] * mm[2] * mm[9] + 
              mm[8] * mm[1] * mm[6] - 
              mm[8] * mm[2] * mm[5];

    float det = mm[0] * inv[0] + mm[1] * inv[4] + mm[2] * inv[8] + mm[3] * inv[12];
    
    if (det == 0.0f) {
        return fmat4_identity(); // Fallback
    }
    
    det = 1.0f / det;
    
    fmat4 result;
    for (int i = 0; i < 16; i++) {
        (&result.matrix.m00)[i] = inv[i] * det;
    }
    
    return result;
}

VECMATH_API dmat2 dmat2_inverse(const dmat2* m)
{
    dmat2 result = { 0 };
    double det = dmat2_determinant(m);

    if (fabs(det) < VECMATH_EPSILON_DZERO) {
        return dmat2_identity();
    }
    
    double inv_det = 1.0f / det;
    result.matrix.m00 =  m->matrix.m11 * inv_det;
    result.matrix.m01 = -m->matrix.m01 * inv_det;
    result.matrix.m10 = -m->matrix.m10 * inv_det;
    result.matrix.m11 =  m->matrix.m00 * inv_det;
    return result;
}

VECMATH_API dmat3 dmat3_inverse(const dmat3* m)
{
    dmat3 result = { 0 };
    double det = dmat3_determinant(m);
    if (fabs(det) < VECMATH_EPSILON_DZERO) {
        return dmat3_identity();
    }
    
    double inv_det = 1.0f / det;
    result.matrix.m00 = (m->matrix.m11 * m->matrix.m22 - m->matrix.m12 * m->matrix.m21) * inv_det;
    result.matrix.m01 = (m->matrix.m02 * m->matrix.m21 - m->matrix.m01 * m->matrix.m22) * inv_det;
    result.matrix.m02 = (m->matrix.m01 * m->matrix.m12 - m->matrix.m02 * m->matrix.m11) * inv_det;
    
    result.matrix.m10 = (m->matrix.m12 * m->matrix.m20 - m->matrix.m10 * m->matrix.m22) * inv_det;
    result.matrix.m11 = (m->matrix.m00 * m->matrix.m22 - m->matrix.m02 * m->matrix.m20) * inv_det;
    result.matrix.m12 = (m->matrix.m02 * m->matrix.m10 - m->matrix.m00 * m->matrix.m12) * inv_det;
    
    result.matrix.m20 = (m->matrix.m10 * m->matrix.m21 - m->matrix.m11 * m->matrix.m20) * inv_det;
    result.matrix.m21 = (m->matrix.m01 * m->matrix.m20 - m->matrix.m00 * m->matrix.m21) * inv_det;
    result.matrix.m22 = (m->matrix.m00 * m->matrix.m11 - m->matrix.m01 * m->matrix.m10) * inv_det;
    return result;
}

VECMATH_API dmat4 dmat4_inverse(const dmat4* m)
{
    dmat4 result = { 0 };
    double inv[4][4];
    
    // calculate cofactors and determinant
    inv[0][0] = m->matrix.m11 * m->matrix.m22 * m->matrix.m33 - 
        m->matrix.m11 * m->matrix.m23 * m->matrix.m32 - 
        m->matrix.m21 * m->matrix.m12 * m->matrix.m33 + 
        m->matrix.m21 * m->matrix.m13 * m->matrix.m32 + 
        m->matrix.m31 * m->matrix.m12 * m->matrix.m23 - 
        m->matrix.m31 * m->matrix.m13 * m->matrix.m22;

    inv[1][0] = -m->matrix.m10 * m->matrix.m22 * m->matrix.m33 + 
        m->matrix.m10 * m->matrix.m23 * m->matrix.m32 + 
        m->matrix.m20 * m->matrix.m12 * m->matrix.m33 - 
        m->matrix.m20 * m->matrix.m13 * m->matrix.m32 - 
        m->matrix.m30 * m->matrix.m12 * m->matrix.m23 + 
        m->matrix.m30 * m->matrix.m13 * m->matrix.m22;

    inv[2][0] = m->matrix.m10 * m->matrix.m21 * m->matrix.m33 - 
        m->matrix.m10 * m->matrix.m23 * m->matrix.m31 - 
        m->matrix.m20 * m->matrix.m11 * m->matrix.m33 + 
        m->matrix.m20 * m->matrix.m13 * m->matrix.m31 + 
        m->matrix.m30 * m->matrix.m11 * m->matrix.m23 - 
        m->matrix.m30 * m->matrix.m13 * m->matrix.m21;

    inv[3][0] = -m->matrix.m10 * m->matrix.m21 * m->matrix.m32 + 
        m->matrix.m10 * m->matrix.m22 * m->matrix.m31 + 
        m->matrix.m20 * m->matrix.m11 * m->matrix.m32 - 
        m->matrix.m20 * m->matrix.m12 * m->matrix.m31 - 
        m->matrix.m30 * m->matrix.m11 * m->matrix.m22 + 
        m->matrix.m30 * m->matrix.m12 * m->matrix.m21;

    inv[0][1] = -m->matrix.m01 * m->matrix.m22 * m->matrix.m33 + 
        m->matrix.m01 * m->matrix.m23 * m->matrix.m32 + 
        m->matrix.m21 * m->matrix.m02 * m->matrix.m33 - 
        m->matrix.m21 * m->matrix.m03 * m->matrix.m32 - 
        m->matrix.m31 * m->matrix.m02 * m->matrix.m23 + 
        m->matrix.m31 * m->matrix.m03 * m->matrix.m22;

    inv[1][1] = m->matrix.m00 * m->matrix.m22 * m->matrix.m33 - 
        m->matrix.m00 * m->matrix.m23 * m->matrix.m32 - 
        m->matrix.m20 * m->matrix.m02 * m->matrix.m33 + 
        m->matrix.m20 * m->matrix.m03 * m->matrix.m32 + 
        m->matrix.m30 * m->matrix.m02 * m->matrix.m23 - 
        m->matrix.m30 * m->matrix.m03 * m->matrix.m22;

    inv[2][1] = -m->matrix.m00 * m->matrix.m21 * m->matrix.m33 + 
        m->matrix.m00 * m->matrix.m23 * m->matrix.m31 + 
        m->matrix.m20 * m->matrix.m01 * m->matrix.m33 - 
        m->matrix.m20 * m->matrix.m03 * m->matrix.m31 - 
        m->matrix.m30 * m->matrix.m01 * m->matrix.m23 + 
        m->matrix.m30 * m->matrix.m03 * m->matrix.m21;

    inv[3][1] = m->matrix.m00 * m->matrix.m21 * m->matrix.m32 - 
        m->matrix.m00 * m->matrix.m22 * m->matrix.m31 - 
        m->matrix.m20 * m->matrix.m01 * m->matrix.m32 + 
        m->matrix.m20 * m->matrix.m02 * m->matrix.m31 + 
        m->matrix.m30 * m->matrix.m01 * m->matrix.m22 - 
        m->matrix.m30 * m->matrix.m02 * m->matrix.m21;

    inv[0][2] = m->matrix.m01 * m->matrix.m12 * m->matrix.m33 - 
        m->matrix.m01 * m->matrix.m13 * m->matrix.m32 - 
        m->matrix.m11 * m->matrix.m02 * m->matrix.m33 + 
        m->matrix.m11 * m->matrix.m03 * m->matrix.m32 + 
        m->matrix.m31 * m->matrix.m02 * m->matrix.m13 - 
        m->matrix.m31 * m->matrix.m03 * m->matrix.m12;

    inv[1][2] = -m->matrix.m00 * m->matrix.m12 * m->matrix.m33 + 
        m->matrix.m00 * m->matrix.m13 * m->matrix.m32 + 
        m->matrix.m10 * m->matrix.m02 * m->matrix.m33 - 
        m->matrix.m10 * m->matrix.m03 * m->matrix.m32 - 
        m->matrix.m30 * m->matrix.m02 * m->matrix.m13 + 
        m->matrix.m30 * m->matrix.m03 * m->matrix.m12;

    inv[2][2] = m->matrix.m00 * m->matrix.m11 * m->matrix.m33 - 
        m->matrix.m00 * m->matrix.m13 * m->matrix.m31 - 
        m->matrix.m10 * m->matrix.m01 * m->matrix.m33 + 
        m->matrix.m10 * m->matrix.m03 * m->matrix.m31 + 
        m->matrix.m30 * m->matrix.m01 * m->matrix.m13 - 
        m->matrix.m30 * m->matrix.m03 * m->matrix.m11;

    inv[3][2] = -m->matrix.m00 * m->matrix.m11 * m->matrix.m32 + 
        m->matrix.m00 * m->matrix.m12 * m->matrix.m31 + 
        m->matrix.m10 * m->matrix.m01 * m->matrix.m32 - 
        m->matrix.m10 * m->matrix.m02 * m->matrix.m31 - 
        m->matrix.m30 * m->matrix.m01 * m->matrix.m12 + 
        m->matrix.m30 * m->matrix.m02 * m->matrix.m11;

    inv[0][3] = -m->matrix.m01 * m->matrix.m12 * m->matrix.m23 + 
        m->matrix.m01 * m->matrix.m13 * m->matrix.m22 + 
        m->matrix.m11 * m->matrix.m02 * m->matrix.m23 - 
        m->matrix.m11 * m->matrix.m03 * m->matrix.m22 - 
        m->matrix.m21 * m->matrix.m02 * m->matrix.m13 + 
        m->matrix.m21 * m->matrix.m03 * m->matrix.m12;

    inv[1][3] = m->matrix.m00 * m->matrix.m12 * m->matrix.m23 - 
        m->matrix.m00 * m->matrix.m13 * m->matrix.m22 - 
        m->matrix.m10 * m->matrix.m02 * m->matrix.m23 + 
        m->matrix.m10 * m->matrix.m03 * m->matrix.m22 + 
        m->matrix.m20 * m->matrix.m02 * m->matrix.m13 - 
        m->matrix.m20 * m->matrix.m03 * m->matrix.m12;

    inv[2][3] = -m->matrix.m00 * m->matrix.m11 * m->matrix.m23 + 
        m->matrix.m00 * m->matrix.m13 * m->matrix.m21 + 
        m->matrix.m10 * m->matrix.m01 * m->matrix.m23 - 
        m->matrix.m10 * m->matrix.m03 * m->matrix.m21 - 
        m->matrix.m20 * m->matrix.m01 * m->matrix.m13 + 
        m->matrix.m20 * m->matrix.m03 * m->matrix.m11;

    inv[3][3] = m->matrix.m00 * m->matrix.m11 * m->matrix.m22 - 
        m->matrix.m00 * m->matrix.m12 * m->matrix.m21 - 
        m->matrix.m10 * m->matrix.m01 * m->matrix.m22 + 
        m->matrix.m10 * m->matrix.m02 * m->matrix.m21 + 
        m->matrix.m20 * m->matrix.m01 * m->matrix.m12 - 
        m->matrix.m20 * m->matrix.m02 * m->matrix.m11;

    // calculate determinant
    double det = m->matrix.m00 * inv[0][0] + m->matrix.m01 * inv[1][0] + m->matrix.m02 * inv[2][0] + m->matrix.m03 * inv[3][0];

    if (fabs(det) < VECMATH_EPSILON_DZERO) {
        return dmat4_identity();
    }

    // scale by 1/determinant
    double inv_det = 1.0f / det;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            result.data[i][j] = inv[i][j] * inv_det;
        }
    }

    return result;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// decompose
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API void fmat4_decompose_rowmajor(const fmat4* m, float3* translation, float3* rotation, float3* scale)
{
    if (translation)    *translation = fmat4_get_translation_rowmajor(m);
    if (scale)          *scale = fmat4_get_scale_rowmajor(m);
    if (rotation)       *rotation = fmat4_get_rotation_rowmajor(m);
}

VECMATH_API void fmat4_decompose_colmajor(const fmat4* m, float3* translation, float3* rotation, float3* scale)
{
    if (translation)    *translation = fmat4_get_translation_colmajor(m);
    if (scale)          *scale = fmat4_get_scale_colmajor(m);
    if (rotation)       *rotation = fmat4_get_rotation_colmajor(m);
}

VECMATH_API void dmat4_decompose_rowmajor(const dmat4* m, double3* translation, double3* rotation, double3* scale)
{
    if (translation)    *translation = dmat4_get_translation_rowmajor(m);
    if (scale)          *scale = dmat4_get_scale_rowmajor(m);
    if (rotation)       *rotation = dmat4_get_rotation_rowmajor(m);
}

VECMATH_API void dmat4_decompose_colmajor(const dmat4* m, double3* translation, double3* rotation, double3* scale)
{
    if (translation)    *translation = dmat4_get_translation_colmajor(m);
    if (scale)          *scale = dmat4_get_scale_colmajor(m);
    if (rotation)       *rotation = dmat4_get_rotation_colmajor(m);
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// get_translation
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API float3 fmat4_get_translation_rowmajor(const fmat4* m)
{
    float3 result = { 0 };
    result.xyz.x = m->matrix.m03;
    result.xyz.y = m->matrix.m13;
    result.xyz.z = m->matrix.m23;
    return result;
}

VECMATH_API float3 fmat4_get_translation_colmajor(const fmat4* m)
{
    float3 result = { 0 };
    result.xyz.x = m->matrix.m30;
    result.xyz.y = m->matrix.m31;
    result.xyz.z = m->matrix.m32;
    return result;
}

VECMATH_API double3 dmat4_get_translation_rowmajor(const dmat4* m)
{
    double3 result = { 0 };
    result.xyz.x = m->matrix.m03;
    result.xyz.y = m->matrix.m13;
    result.xyz.z = m->matrix.m23;
    return result;
}

VECMATH_API double3 dmat4_get_translation_colmajor(const dmat4* m)
{
    double3 result = { 0 };
    result.xyz.x = m->matrix.m30;
    result.xyz.y = m->matrix.m31;
    result.xyz.z = m->matrix.m32;
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// get_rotation
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API float3 fmat4_get_rotation_rowmajor(const fmat4* m)
{
    float3 scale = fmat4_get_scale_rowmajor(m);
    float3 rotation = { 0 };
    
    if (scale.xyz.x < VECMATH_EPSILON_FZERO || scale.xyz.y < VECMATH_EPSILON_FZERO || scale.xyz.z < VECMATH_EPSILON_FZERO) {
        return rotation;
    }
    
    // remove scale from rows
    float m00 = m->matrix.m00 / scale.xyz.x;
    float m01 = m->matrix.m01 / scale.xyz.x;
    float m02 = m->matrix.m02 / scale.xyz.x;
    
    float m10 = m->matrix.m10 / scale.xyz.y;
    float m11 = m->matrix.m11 / scale.xyz.y;
    float m12 = m->matrix.m12 / scale.xyz.y;
    
    float m20 = m->matrix.m20 / scale.xyz.z;
    float m21 = m->matrix.m21 / scale.xyz.z;
    float m22 = m->matrix.m22 / scale.xyz.z;
    
    // extract Euler angles from row-major rotation matrix
    rotation.xyz.y = atan2f(m02, sqrtf(m00*m00 + m01*m01));
    
    if (fabsf(rotation.xyz.y - VECMATH_EPSILON_PI / 2.0f) < VECMATH_EPSILON_FZERO) {
        rotation.xyz.x = atan2f(m10, m11);
        rotation.xyz.z = 0.0f;
    }

    else if (fabsf(rotation.xyz.y + VECMATH_EPSILON_PI / 2.0f) < VECMATH_EPSILON_FZERO) {
        rotation.xyz.x = -atan2f(m10, m11);
        rotation.xyz.z = 0.0f;
    }

    else {
        rotation.xyz.x = atan2f(m12, m22);
        rotation.xyz.z = atan2f(m01, m00);
    }
    
    return rotation;
}

VECMATH_API float3 fmat4_get_rotation_colmajor(const fmat4* m)
{
    float3 scale = fmat4_get_scale_colmajor(m);
    float3 rotation = { 0 };
    
    if (scale.xyz.x < VECMATH_EPSILON_FZERO || scale.xyz.y < VECMATH_EPSILON_FZERO || scale.xyz.z < VECMATH_EPSILON_FZERO) {
        return rotation;
    }
    
    // remove scale from columns
    float m00 = m->matrix.m00 / scale.xyz.x;
    float m10 = m->matrix.m10 / scale.xyz.x;
    float m20 = m->matrix.m20 / scale.xyz.x;
    
    float m01 = m->matrix.m01 / scale.xyz.y;
    float m11 = m->matrix.m11 / scale.xyz.y;
    float m21 = m->matrix.m21 / scale.xyz.y;
    
    float m02 = m->matrix.m02 / scale.xyz.z;
    float m12 = m->matrix.m12 / scale.xyz.z;
    float m22 = m->matrix.m22 / scale.xyz.z;
    
    // extract euler angles from column-major rotation matrix
    rotation.xyz.y = atan2f(-m20, sqrtf(m00*m00 + m10*m10));
    
    if (fabsf(rotation.xyz.y - VECMATH_EPSILON_PI / 2.0f) < VECMATH_EPSILON_FZERO) {
        rotation.xyz.x = atan2f(m01, m11);
        rotation.xyz.z = 0.0f;
    }

    else if (fabsf(rotation.xyz.y + VECMATH_EPSILON_PI / 2.0f) < VECMATH_EPSILON_FZERO) {
        rotation.xyz.x = -atan2f(m01, m11);
        rotation.xyz.z = 0.0f;
    }

    else {
        rotation.xyz.x = atan2f(m21, m22);
        rotation.xyz.z = atan2f(m10, m00);
    }
    
    return rotation;
}

VECMATH_API double3 dmat4_get_rotation_rowmajor(const dmat4* m)
{
    double3 scale = dmat4_get_scale_rowmajor(m);
    double3 rotation = { 0 };
    
    if (scale.xyz.x < VECMATH_EPSILON_DZERO || scale.xyz.y < VECMATH_EPSILON_DZERO || scale.xyz.z < VECMATH_EPSILON_DZERO) {
        return rotation;
    }
    
    // remove scale from rows
    double m00 = m->matrix.m00 / scale.xyz.x;
    double m01 = m->matrix.m01 / scale.xyz.x;
    double m02 = m->matrix.m02 / scale.xyz.x;
    
    double m10 = m->matrix.m10 / scale.xyz.y;
    double m11 = m->matrix.m11 / scale.xyz.y;
    double m12 = m->matrix.m12 / scale.xyz.y;
    
    double m20 = m->matrix.m20 / scale.xyz.z;
    double m21 = m->matrix.m21 / scale.xyz.z;
    double m22 = m->matrix.m22 / scale.xyz.z;
    
    // extract Euler angles from row-major rotation matrix
    rotation.xyz.y = atan2(m02, sqrt(m00*m00 + m01*m01));
    
    if (fabs(rotation.xyz.y - VECMATH_EPSILON_PI / 2.0) < VECMATH_EPSILON_DZERO) {
        rotation.xyz.x = atan2(m10, m11);
        rotation.xyz.z = 0.0;
    }

    else if (fabs(rotation.xyz.y + VECMATH_EPSILON_PI / 2.0) < VECMATH_EPSILON_DZERO) {
        rotation.xyz.x = -atan2(m10, m11);
        rotation.xyz.z = 0.0;
    }

    else {
        rotation.xyz.x = atan2(m12, m22);
        rotation.xyz.z = atan2(m01, m00);
    }
    
    return rotation;
}

VECMATH_API double3 dmat4_get_rotation_colmajor(const dmat4* m)
{
    double3 scale = dmat4_get_scale_colmajor(m);
    double3 rotation = { 0 };
    
    if (scale.xyz.x < VECMATH_EPSILON_DZERO || scale.xyz.y < VECMATH_EPSILON_DZERO || scale.xyz.z < VECMATH_EPSILON_DZERO) {
        return rotation;
    }
    
    // remove scale from columns
    double m00 = m->matrix.m00 / scale.xyz.x;
    double m10 = m->matrix.m10 / scale.xyz.x;
    double m20 = m->matrix.m20 / scale.xyz.x;
    
    double m01 = m->matrix.m01 / scale.xyz.y;
    double m11 = m->matrix.m11 / scale.xyz.y;
    double m21 = m->matrix.m21 / scale.xyz.y;
    
    double m02 = m->matrix.m02 / scale.xyz.z;
    double m12 = m->matrix.m12 / scale.xyz.z;
    double m22 = m->matrix.m22 / scale.xyz.z;
    
    // extract euler angles from column-major rotation matrix
    rotation.xyz.y = atan2(-m20, sqrt(m00*m00 + m10*m10));
    
    if (fabs(rotation.xyz.y - VECMATH_EPSILON_PI / 2.0) < VECMATH_EPSILON_DZERO) {
        rotation.xyz.x = atan2(m01, m11);
        rotation.xyz.z = 0.0;
    }

    else if (fabs(rotation.xyz.y + VECMATH_EPSILON_PI / 2.0) < VECMATH_EPSILON_DZERO) {
        rotation.xyz.x = -atan2(m01, m11);
        rotation.xyz.z = 0.0;
    }

    else {
        rotation.xyz.x = atan2(m21, m22);
        rotation.xyz.z = atan2(m10, m00);
    }
    
    return rotation;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// get_scale
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API float3 fmat4_get_scale_rowmajor(const fmat4* m)
{
    // scale is length of each row (basis vectors)
    float3 scale;
    scale.xyz.x = sqrtf(m->matrix.m00 * m->matrix.m00 + m->matrix.m01 * m->matrix.m01 + m->matrix.m02*m->matrix.m02);
    scale.xyz.y = sqrtf(m->matrix.m10 * m->matrix.m10 + m->matrix.m11 * m->matrix.m11 + m->matrix.m12*m->matrix.m12);
    scale.xyz.z = sqrtf(m->matrix.m20 * m->matrix.m20 + m->matrix.m21 * m->matrix.m21 + m->matrix.m22*m->matrix.m22);
    
    // check for negative scale
    float det = 
        m->matrix.m00 * (m->matrix.m11 * m->matrix.m22 - m->matrix.m12 * m->matrix.m21) -
        m->matrix.m01 * (m->matrix.m10 * m->matrix.m22 - m->matrix.m12 * m->matrix.m20) +
        m->matrix.m02 * (m->matrix.m10 * m->matrix.m21 - m->matrix.m11 * m->matrix.m20);
    
    if (det < 0) {
        scale.xyz.x = -scale.xyz.x;
    }
    
    return scale;
}

VECMATH_API float3 fmat4_get_scale_colmajor(const fmat4* m)
{
    // scale is length of each column (basis vectors)
    float3 scale;
    scale.xyz.x = sqrtf(m->matrix.m00 * m->matrix.m00 + m->matrix.m10 * m->matrix.m10 + m->matrix.m20 * m->matrix.m20);
    scale.xyz.y = sqrtf(m->matrix.m01 * m->matrix.m01 + m->matrix.m11 * m->matrix.m11 + m->matrix.m21 * m->matrix.m21);
    scale.xyz.z = sqrtf(m->matrix.m02 * m->matrix.m02 + m->matrix.m12 * m->matrix.m12 + m->matrix.m22 * m->matrix.m22);
    
    // check for negative scale
    float det =
        m->matrix.m00 * (m->matrix.m11 * m->matrix.m22 - m->matrix.m12 * m->matrix.m21) -
        m->matrix.m01 * (m->matrix.m10 * m->matrix.m22 - m->matrix.m12 * m->matrix.m20) +
        m->matrix.m02 * (m->matrix.m10 * m->matrix.m21 - m->matrix.m11 * m->matrix.m20);
    
    if (det < 0) {
        scale.xyz.x = -scale.xyz.x;
    }
    
    return scale;
}

VECMATH_API double3 dmat4_get_scale_rowmajor(const dmat4* m)
{
    // scale is length of each row (basis vectors)
    double3 scale;
    scale.xyz.x = sqrtf(m->matrix.m00 * m->matrix.m00 + m->matrix.m01 * m->matrix.m01 + m->matrix.m02*m->matrix.m02);
    scale.xyz.y = sqrtf(m->matrix.m10 * m->matrix.m10 + m->matrix.m11 * m->matrix.m11 + m->matrix.m12*m->matrix.m12);
    scale.xyz.z = sqrtf(m->matrix.m20 * m->matrix.m20 + m->matrix.m21 * m->matrix.m21 + m->matrix.m22*m->matrix.m22);
    
    // check for negative scale
    double det = 
        m->matrix.m00 * (m->matrix.m11 * m->matrix.m22 - m->matrix.m12 * m->matrix.m21) -
        m->matrix.m01 * (m->matrix.m10 * m->matrix.m22 - m->matrix.m12 * m->matrix.m20) +
        m->matrix.m02 * (m->matrix.m10 * m->matrix.m21 - m->matrix.m11 * m->matrix.m20);
    
    if (det < 0) {
        scale.xyz.x = -scale.xyz.x;
    }
    
    return scale;
}

VECMATH_API double3 dmat4_get_scale_colmajor(const dmat4* m)
{
    // scale is length of each column (basis vectors)
    double3 scale;
    scale.xyz.x = sqrtf(m->matrix.m00 * m->matrix.m00 + m->matrix.m10 * m->matrix.m10 + m->matrix.m20 * m->matrix.m20);
    scale.xyz.y = sqrtf(m->matrix.m01 * m->matrix.m01 + m->matrix.m11 * m->matrix.m11 + m->matrix.m21 * m->matrix.m21);
    scale.xyz.z = sqrtf(m->matrix.m02 * m->matrix.m02 + m->matrix.m12 * m->matrix.m12 + m->matrix.m22 * m->matrix.m22);
    
    // check for negative scale
    double det =
        m->matrix.m00 * (m->matrix.m11 * m->matrix.m22 - m->matrix.m12 * m->matrix.m21) -
        m->matrix.m01 * (m->matrix.m10 * m->matrix.m22 - m->matrix.m12 * m->matrix.m20) +
        m->matrix.m02 * (m->matrix.m10 * m->matrix.m21 - m->matrix.m11 * m->matrix.m20);
    
    if (det < 0) {
        scale.xyz.x = -scale.xyz.x;
    }
    
    return scale;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// translate
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API fmat4 fmat4_translate_rowmajor(const fmat4* m, const float3* dir)
{
    fmat4 result = *m;

    result.data[3][0] += dir->xyz.x;
    result.data[3][1] += dir->xyz.y;
    result.data[3][2] += dir->xyz.z;
    return result;
}

VECMATH_API fmat4 fmat4_translate_colmajor(const fmat4* m, const float3 *dir)
{
    fmat4 result = *m;

    // column-major: translation goes in m03, m13, m23
    result.matrix.m03 = m->matrix.m00 * dir->xyz.x + m->matrix.m01 * dir->xyz.y + m->matrix.m02 * dir->xyz.z + m->matrix.m03;
    result.matrix.m13 = m->matrix.m10 * dir->xyz.x + m->matrix.m11 * dir->xyz.y + m->matrix.m12 * dir->xyz.z + m->matrix.m13;
    result.matrix.m23 = m->matrix.m20 * dir->xyz.x + m->matrix.m21 * dir->xyz.y + m->matrix.m22 * dir->xyz.z + m->matrix.m23;
    result.matrix.m33 = m->matrix.m30 * dir->xyz.x + m->matrix.m31 * dir->xyz.y + m->matrix.m32 * dir->xyz.z + m->matrix.m33;

    return result;
}

VECMATH_API dmat4 dmat4_translate_rowmajor(const dmat4* m, const double3* dir)
{
    dmat4 result = *m;

    result.data[3][0] += dir->xyz.x;
    result.data[3][1] += dir->xyz.y;
    result.data[3][2] += dir->xyz.z;
    return result;
}

VECMATH_API dmat4 dmat4_translate_colmajor(const dmat4* m, const double3* dir)
{
    dmat4 result = *m;

    result.data[0][3] += dir->xyz.x;
    result.data[1][3] += dir->xyz.y;
    result.data[2][3] += dir->xyz.z;
    return result;
}


/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// rotate
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API fmat4 fmat4_rotate_colmajor(const fmat4* m, float angle, const float3* axis)
{
    // normalize axis and compute temp values
    float3 axis_n = float3_normalize(axis);
    float c = cosf(angle);
    float s = sinf(angle);
    float one_minus_c = 1.0f - c;
    fmat4 rotate = { 0 };
    
    // column 0
    rotate.data[0][0] = c + axis_n.xyz.x * axis_n.xyz.x * one_minus_c;
    rotate.data[1][0] = axis_n.xyz.y * axis_n.xyz.x * one_minus_c + axis_n.xyz.z * s;
    rotate.data[2][0] = axis_n.xyz.z * axis_n.xyz.x * one_minus_c - axis_n.xyz.y * s;
    // column 1
    rotate.data[0][1] = axis_n.xyz.x * axis_n.xyz.y * one_minus_c - axis_n.xyz.z * s;
    rotate.data[1][1] = c + axis_n.xyz.y * axis_n.xyz.y * one_minus_c;
    rotate.data[2][1] = axis_n.xyz.z * axis_n.xyz.y * one_minus_c + axis_n.xyz.x * s;
    // column 2
    rotate.data[0][2] = axis_n.xyz.x * axis_n.xyz.z * one_minus_c + axis_n.xyz.y * s;
    rotate.data[1][2] = axis_n.xyz.y * axis_n.xyz.z * one_minus_c - axis_n.xyz.x * s;
    rotate.data[2][2] = c + axis_n.xyz.z * axis_n.xyz.z * one_minus_c;
    
    // column 3 (identity)
    rotate.data[3][3] = 1.0f;
    
    // multiply: result = m * rotate (column-major)
    fmat4 result = { 0 };
    
    // column 0
    result.data[0][0] = m->data[0][0]*rotate.data[0][0] + m->data[0][1]*rotate.data[1][0] + m->data[0][2]*rotate.data[2][0] + m->data[0][3]*rotate.data[3][0];
    result.data[1][0] = m->data[1][0]*rotate.data[0][0] + m->data[1][1]*rotate.data[1][0] + m->data[1][2]*rotate.data[2][0] + m->data[1][3]*rotate.data[3][0];
    result.data[2][0] = m->data[2][0]*rotate.data[0][0] + m->data[2][1]*rotate.data[1][0] + m->data[2][2]*rotate.data[2][0] + m->data[2][3]*rotate.data[3][0];
    result.data[3][0] = m->data[3][0]*rotate.data[0][0] + m->data[3][1]*rotate.data[1][0] + m->data[3][2]*rotate.data[2][0] + m->data[3][3]*rotate.data[3][0];
    // column 1
    result.data[0][1] = m->data[0][0]*rotate.data[0][1] + m->data[0][1]*rotate.data[1][1] + m->data[0][2]*rotate.data[2][1] + m->data[0][3]*rotate.data[3][1];
    result.data[1][1] = m->data[1][0]*rotate.data[0][1] + m->data[1][1]*rotate.data[1][1] + m->data[1][2]*rotate.data[2][1] + m->data[1][3]*rotate.data[3][1];
    result.data[2][1] = m->data[2][0]*rotate.data[0][1] + m->data[2][1]*rotate.data[1][1] + m->data[2][2]*rotate.data[2][1] + m->data[2][3]*rotate.data[3][1];
    result.data[3][1] = m->data[3][0]*rotate.data[0][1] + m->data[3][1]*rotate.data[1][1] + m->data[3][2]*rotate.data[2][1] + m->data[3][3]*rotate.data[3][1];
    // column 2
    result.data[0][2] = m->data[0][0]*rotate.data[0][2] + m->data[0][1]*rotate.data[1][2] + m->data[0][2]*rotate.data[2][2] + m->data[0][3]*rotate.data[3][2];
    result.data[1][2] = m->data[1][0]*rotate.data[0][2] + m->data[1][1]*rotate.data[1][2] + m->data[1][2]*rotate.data[2][2] + m->data[1][3]*rotate.data[3][2];
    result.data[2][2] = m->data[2][0]*rotate.data[0][2] + m->data[2][1]*rotate.data[1][2] + m->data[2][2]*rotate.data[2][2] + m->data[2][3]*rotate.data[3][2];
    result.data[3][2] = m->data[3][0]*rotate.data[0][2] + m->data[3][1]*rotate.data[1][2] + m->data[3][2]*rotate.data[2][2] + m->data[3][3]*rotate.data[3][2];
    // column 3
    result.data[0][3] = m->data[0][3];
    result.data[1][3] = m->data[1][3];
    result.data[2][3] = m->data[2][3];
    result.data[3][3] = m->data[3][3];
    
    return result;
}

VECMATH_API fmat4 fmat4_rotate_rowmajor(const fmat4* m, float angle, const float3* axis)
{
    // normalize axis and compute temp values
    float3 axis_n = float3_normalize(axis);
    float c = cosf(angle);
    float s = sinf(angle);
    float one_minus_c = 1.0f - c;
    
    // construct rotation matrix (row-major)
    fmat4 rotate = { 0 };
    
    // row 0
    rotate.data[0][0] = c + axis_n.xyz.x * axis_n.xyz.x * one_minus_c;
    rotate.data[0][1] = axis_n.xyz.x * axis_n.xyz.y * one_minus_c - axis_n.xyz.z * s;
    rotate.data[0][2] = axis_n.xyz.x * axis_n.xyz.z * one_minus_c + axis_n.xyz.y * s;
    rotate.data[0][3] = 0.0f;
    // row 1
    rotate.data[1][0] = axis_n.xyz.y * axis_n.xyz.x * one_minus_c + axis_n.xyz.z * s;
    rotate.data[1][1] = c + axis_n.xyz.y * axis_n.xyz.y * one_minus_c;
    rotate.data[1][2] = axis_n.xyz.y * axis_n.xyz.z * one_minus_c - axis_n.xyz.x * s;
    rotate.data[1][3] = 0.0f;
    // row 2
    rotate.data[2][0] = axis_n.xyz.z * axis_n.xyz.x * one_minus_c - axis_n.xyz.y * s;
    rotate.data[2][1] = axis_n.xyz.z * axis_n.xyz.y * one_minus_c + axis_n.xyz.x * s;
    rotate.data[2][2] = c + axis_n.xyz.z * axis_n.xyz.z * one_minus_c;
    rotate.data[2][3] = 0.0f;
    // row 3
    rotate.data[3][0] = 0.0f;
    rotate.data[3][1] = 0.0f;
    rotate.data[3][2] = 0.0f;
    rotate.data[3][3] = 1.0f;
    
    // multiply: result = rotate * m (row-major)
    fmat4 result = { 0 };
    
    // row 0
    result.data[0][0] = rotate.data[0][0]*m->data[0][0] + rotate.data[0][1]*m->data[1][0] + rotate.data[0][2]*m->data[2][0] + rotate.data[0][3]*m->data[3][0];
    result.data[0][1] = rotate.data[0][0]*m->data[0][1] + rotate.data[0][1]*m->data[1][1] + rotate.data[0][2]*m->data[2][1] + rotate.data[0][3]*m->data[3][1];
    result.data[0][2] = rotate.data[0][0]*m->data[0][2] + rotate.data[0][1]*m->data[1][2] + rotate.data[0][2]*m->data[2][2] + rotate.data[0][3]*m->data[3][2];
    result.data[0][3] = rotate.data[0][0]*m->data[0][3] + rotate.data[0][1]*m->data[1][3] + rotate.data[0][2]*m->data[2][3] + rotate.data[0][3]*m->data[3][3];
    // row 1
    result.data[1][0] = rotate.data[1][0]*m->data[0][0] + rotate.data[1][1]*m->data[1][0] + rotate.data[1][2]*m->data[2][0] + rotate.data[1][3]*m->data[3][0];
    result.data[1][1] = rotate.data[1][0]*m->data[0][1] + rotate.data[1][1]*m->data[1][1] + rotate.data[1][2]*m->data[2][1] + rotate.data[1][3]*m->data[3][1];
    result.data[1][2] = rotate.data[1][0]*m->data[0][2] + rotate.data[1][1]*m->data[1][2] + rotate.data[1][2]*m->data[2][2] + rotate.data[1][3]*m->data[3][2];
    result.data[1][3] = rotate.data[1][0]*m->data[0][3] + rotate.data[1][1]*m->data[1][3] + rotate.data[1][2]*m->data[2][3] + rotate.data[1][3]*m->data[3][3];
    // row 2
    result.data[2][0] = rotate.data[2][0]*m->data[0][0] + rotate.data[2][1]*m->data[1][0] + rotate.data[2][2]*m->data[2][0] + rotate.data[2][3]*m->data[3][0];
    result.data[2][1] = rotate.data[2][0]*m->data[0][1] + rotate.data[2][1]*m->data[1][1] + rotate.data[2][2]*m->data[2][1] + rotate.data[2][3]*m->data[3][1];
    result.data[2][2] = rotate.data[2][0]*m->data[0][2] + rotate.data[2][1]*m->data[1][2] + rotate.data[2][2]*m->data[2][2] + rotate.data[2][3]*m->data[3][2];
    result.data[2][3] = rotate.data[2][0]*m->data[0][3] + rotate.data[2][1]*m->data[1][3] + rotate.data[2][2]*m->data[2][3] + rotate.data[2][3]*m->data[3][3];
    // row 3 
    result.data[3][0] = m->data[3][0];
    result.data[3][1] = m->data[3][1];
    result.data[3][2] = m->data[3][2];
    result.data[3][3] = m->data[3][3];
    
    return result;
}

VECMATH_API dmat4 dmat4_rotate_colmajor(const dmat4* m, double angle, const double3* axis)
{
    // normalize axis and compute temp values
    double3 axis_n = double3_normalize(axis);
    double c = cos(angle);
    double s = sin(angle);
    double one_minus_c = 1.0 - c;
    dmat4 rotate = { 0 };
    
    // column 0
    rotate.data[0][0] = c + axis_n.xyz.x * axis_n.xyz.x * one_minus_c;
    rotate.data[1][0] = axis_n.xyz.y * axis_n.xyz.x * one_minus_c + axis_n.xyz.z * s;
    rotate.data[2][0] = axis_n.xyz.z * axis_n.xyz.x * one_minus_c - axis_n.xyz.y * s;
    // column 1
    rotate.data[0][1] = axis_n.xyz.x * axis_n.xyz.y * one_minus_c - axis_n.xyz.z * s;
    rotate.data[1][1] = c + axis_n.xyz.y * axis_n.xyz.y * one_minus_c;
    rotate.data[2][1] = axis_n.xyz.z * axis_n.xyz.y * one_minus_c + axis_n.xyz.x * s;
    // column 2
    rotate.data[0][2] = axis_n.xyz.x * axis_n.xyz.z * one_minus_c + axis_n.xyz.y * s;
    rotate.data[1][2] = axis_n.xyz.y * axis_n.xyz.z * one_minus_c - axis_n.xyz.x * s;
    rotate.data[2][2] = c + axis_n.xyz.z * axis_n.xyz.z * one_minus_c;
    
    // column 3 (identity)
    rotate.data[3][3] = 1.0f;
    
    // multiply: result = m * rotate (column-major)
    dmat4 result = { 0 };
    
    // column 0
    result.data[0][0] = m->data[0][0]*rotate.data[0][0] + m->data[0][1]*rotate.data[1][0] + m->data[0][2]*rotate.data[2][0] + m->data[0][3]*rotate.data[3][0];
    result.data[1][0] = m->data[1][0]*rotate.data[0][0] + m->data[1][1]*rotate.data[1][0] + m->data[1][2]*rotate.data[2][0] + m->data[1][3]*rotate.data[3][0];
    result.data[2][0] = m->data[2][0]*rotate.data[0][0] + m->data[2][1]*rotate.data[1][0] + m->data[2][2]*rotate.data[2][0] + m->data[2][3]*rotate.data[3][0];
    result.data[3][0] = m->data[3][0]*rotate.data[0][0] + m->data[3][1]*rotate.data[1][0] + m->data[3][2]*rotate.data[2][0] + m->data[3][3]*rotate.data[3][0];
    // column 1
    result.data[0][1] = m->data[0][0]*rotate.data[0][1] + m->data[0][1]*rotate.data[1][1] + m->data[0][2]*rotate.data[2][1] + m->data[0][3]*rotate.data[3][1];
    result.data[1][1] = m->data[1][0]*rotate.data[0][1] + m->data[1][1]*rotate.data[1][1] + m->data[1][2]*rotate.data[2][1] + m->data[1][3]*rotate.data[3][1];
    result.data[2][1] = m->data[2][0]*rotate.data[0][1] + m->data[2][1]*rotate.data[1][1] + m->data[2][2]*rotate.data[2][1] + m->data[2][3]*rotate.data[3][1];
    result.data[3][1] = m->data[3][0]*rotate.data[0][1] + m->data[3][1]*rotate.data[1][1] + m->data[3][2]*rotate.data[2][1] + m->data[3][3]*rotate.data[3][1];
    // column 2
    result.data[0][2] = m->data[0][0]*rotate.data[0][2] + m->data[0][1]*rotate.data[1][2] + m->data[0][2]*rotate.data[2][2] + m->data[0][3]*rotate.data[3][2];
    result.data[1][2] = m->data[1][0]*rotate.data[0][2] + m->data[1][1]*rotate.data[1][2] + m->data[1][2]*rotate.data[2][2] + m->data[1][3]*rotate.data[3][2];
    result.data[2][2] = m->data[2][0]*rotate.data[0][2] + m->data[2][1]*rotate.data[1][2] + m->data[2][2]*rotate.data[2][2] + m->data[2][3]*rotate.data[3][2];
    result.data[3][2] = m->data[3][0]*rotate.data[0][2] + m->data[3][1]*rotate.data[1][2] + m->data[3][2]*rotate.data[2][2] + m->data[3][3]*rotate.data[3][2];
    // column 3
    result.data[0][3] = m->data[0][3];
    result.data[1][3] = m->data[1][3];
    result.data[2][3] = m->data[2][3];
    result.data[3][3] = m->data[3][3];
    
    return result;
}

VECMATH_API dmat4 dmat4_rotate_rowmajor(const dmat4* m, double angle, const double3* axis)
{
    // normalize axis and compute temp values
    double3 axis_n = double3_normalize(axis);
    double c = cos(angle);
    double s = sin(angle);
    double one_minus_c = 1.0f - c;
    
    // construct rotation matrix (row-major)
    fmat4 rotate = { 0 };
    
    // row 0
    rotate.data[0][0] = c + axis_n.xyz.x * axis_n.xyz.x * one_minus_c;
    rotate.data[0][1] = axis_n.xyz.x * axis_n.xyz.y * one_minus_c - axis_n.xyz.z * s;
    rotate.data[0][2] = axis_n.xyz.x * axis_n.xyz.z * one_minus_c + axis_n.xyz.y * s;
    rotate.data[0][3] = 0.0f;
    // row 1
    rotate.data[1][0] = axis_n.xyz.y * axis_n.xyz.x * one_minus_c + axis_n.xyz.z * s;
    rotate.data[1][1] = c + axis_n.xyz.y * axis_n.xyz.y * one_minus_c;
    rotate.data[1][2] = axis_n.xyz.y * axis_n.xyz.z * one_minus_c - axis_n.xyz.x * s;
    rotate.data[1][3] = 0.0f;
    // row 2
    rotate.data[2][0] = axis_n.xyz.z * axis_n.xyz.x * one_minus_c - axis_n.xyz.y * s;
    rotate.data[2][1] = axis_n.xyz.z * axis_n.xyz.y * one_minus_c + axis_n.xyz.x * s;
    rotate.data[2][2] = c + axis_n.xyz.z * axis_n.xyz.z * one_minus_c;
    rotate.data[2][3] = 0.0f;
    // row 3
    rotate.data[3][0] = 0.0f;
    rotate.data[3][1] = 0.0f;
    rotate.data[3][2] = 0.0f;
    rotate.data[3][3] = 1.0f;
    
    // multiply: result = rotate * m (row-major)
    dmat4 result = { 0 };
    
    // row 0
    result.data[0][0] = rotate.data[0][0]*m->data[0][0] + rotate.data[0][1]*m->data[1][0] + rotate.data[0][2]*m->data[2][0] + rotate.data[0][3]*m->data[3][0];
    result.data[0][1] = rotate.data[0][0]*m->data[0][1] + rotate.data[0][1]*m->data[1][1] + rotate.data[0][2]*m->data[2][1] + rotate.data[0][3]*m->data[3][1];
    result.data[0][2] = rotate.data[0][0]*m->data[0][2] + rotate.data[0][1]*m->data[1][2] + rotate.data[0][2]*m->data[2][2] + rotate.data[0][3]*m->data[3][2];
    result.data[0][3] = rotate.data[0][0]*m->data[0][3] + rotate.data[0][1]*m->data[1][3] + rotate.data[0][2]*m->data[2][3] + rotate.data[0][3]*m->data[3][3];
    // row 1
    result.data[1][0] = rotate.data[1][0]*m->data[0][0] + rotate.data[1][1]*m->data[1][0] + rotate.data[1][2]*m->data[2][0] + rotate.data[1][3]*m->data[3][0];
    result.data[1][1] = rotate.data[1][0]*m->data[0][1] + rotate.data[1][1]*m->data[1][1] + rotate.data[1][2]*m->data[2][1] + rotate.data[1][3]*m->data[3][1];
    result.data[1][2] = rotate.data[1][0]*m->data[0][2] + rotate.data[1][1]*m->data[1][2] + rotate.data[1][2]*m->data[2][2] + rotate.data[1][3]*m->data[3][2];
    result.data[1][3] = rotate.data[1][0]*m->data[0][3] + rotate.data[1][1]*m->data[1][3] + rotate.data[1][2]*m->data[2][3] + rotate.data[1][3]*m->data[3][3];
    // row 2
    result.data[2][0] = rotate.data[2][0]*m->data[0][0] + rotate.data[2][1]*m->data[1][0] + rotate.data[2][2]*m->data[2][0] + rotate.data[2][3]*m->data[3][0];
    result.data[2][1] = rotate.data[2][0]*m->data[0][1] + rotate.data[2][1]*m->data[1][1] + rotate.data[2][2]*m->data[2][1] + rotate.data[2][3]*m->data[3][1];
    result.data[2][2] = rotate.data[2][0]*m->data[0][2] + rotate.data[2][1]*m->data[1][2] + rotate.data[2][2]*m->data[2][2] + rotate.data[2][3]*m->data[3][2];
    result.data[2][3] = rotate.data[2][0]*m->data[0][3] + rotate.data[2][1]*m->data[1][3] + rotate.data[2][2]*m->data[2][3] + rotate.data[2][3]*m->data[3][3];
    // row 3 
    result.data[3][0] = m->data[3][0];
    result.data[3][1] = m->data[3][1];
    result.data[3][2] = m->data[3][2];
    result.data[3][3] = m->data[3][3];
    
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// scale
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API fmat4 fmat4_scale_rowmajor(const fmat4* m, const float3* dim)
{
    fmat4 result = { 0 };

    result.data[0][0] = m->data[0][0] * dim->xyz.x;
    result.data[0][1] = m->data[0][1] * dim->xyz.x;
    result.data[0][2] = m->data[0][2] * dim->xyz.x;
    result.data[0][3] = m->data[0][3] * dim->xyz.x;

    result.data[1][0] = m->data[1][0] * dim->xyz.y;
    result.data[1][1] = m->data[1][1] * dim->xyz.y;
    result.data[1][2] = m->data[1][2] * dim->xyz.y;
    result.data[1][3] = m->data[1][3] * dim->xyz.y;

    result.data[2][0] = m->data[2][0] * dim->xyz.z;
    result.data[2][1] = m->data[2][1] * dim->xyz.z;
    result.data[2][2] = m->data[2][2] * dim->xyz.z;
    result.data[2][3] = m->data[2][3] * dim->xyz.z;

    result.data[3][0] = m->data[3][0];
    result.data[3][1] = m->data[3][1];
    result.data[3][2] = m->data[3][2];
    result.data[3][3] = m->data[3][3];
    return result;
}

VECMATH_API fmat4 fmat4_scale_colmajor(const fmat4* m, const float3* dim)
{
    fmat4 result = *m;

    // column 0 (X axis)
    result.matrix.m00 *= dim->xyz.x;
    result.matrix.m10 *= dim->xyz.x;
    result.matrix.m20 *= dim->xyz.x;
    result.matrix.m30 *= dim->xyz.x;

    // column 1 (Y axis)
    result.matrix.m01 *= dim->xyz.y;
    result.matrix.m11 *= dim->xyz.y;
    result.matrix.m21 *= dim->xyz.y;
    result.matrix.m31 *= dim->xyz.y;

    // column 2 (Z axis)
    result.matrix.m02 *= dim->xyz.z;
    result.matrix.m12 *= dim->xyz.z;
    result.matrix.m22 *= dim->xyz.z;
    result.matrix.m32 *= dim->xyz.z;

    return result;
}

VECMATH_API dmat4 dmat4_scale_rowmajor(const dmat4* m, const double3* dim)
{
    dmat4 result = { 0 };

    result.data[0][0] = m->data[0][0] * dim->xyz.x;
    result.data[0][1] = m->data[0][1] * dim->xyz.x;
    result.data[0][2] = m->data[0][2] * dim->xyz.x;
    result.data[0][3] = m->data[0][3] * dim->xyz.x;

    result.data[1][0] = m->data[1][0] * dim->xyz.y;
    result.data[1][1] = m->data[1][1] * dim->xyz.y;
    result.data[1][2] = m->data[1][2] * dim->xyz.y;
    result.data[1][3] = m->data[1][3] * dim->xyz.y;

    result.data[2][0] = m->data[2][0] * dim->xyz.z;
    result.data[2][1] = m->data[2][1] * dim->xyz.z;
    result.data[2][2] = m->data[2][2] * dim->xyz.z;
    result.data[2][3] = m->data[2][3] * dim->xyz.z;

    result.data[3][0] = m->data[3][0];
    result.data[3][1] = m->data[3][1];
    result.data[3][2] = m->data[3][2];
    result.data[3][3] = m->data[3][3];
    return result;
}

VECMATH_API dmat4 dmat4_scale_colmajor(const dmat4* m, const double3* dim)
{
    dmat4 result = { 0 };

    result.data[0][0] = m->data[0][0] * dim->xyz.x;
    result.data[1][0] = m->data[1][0] * dim->xyz.x;
    result.data[2][0] = m->data[2][0] * dim->xyz.x;
    result.data[3][0] = m->data[3][0] * dim->xyz.x;

    result.data[0][1] = m->data[0][1] * dim->xyz.y;
    result.data[1][1] = m->data[1][1] * dim->xyz.y;
    result.data[2][1] = m->data[2][1] * dim->xyz.y;
    result.data[3][1] = m->data[3][1] * dim->xyz.y;

    result.data[0][2] = m->data[0][2] * dim->xyz.z;
    result.data[1][2] = m->data[1][2] * dim->xyz.z;
    result.data[2][2] = m->data[2][2] * dim->xyz.z;
    result.data[3][2] = m->data[3][2] * dim->xyz.z;

    result.data[0][3] = m->data[0][3];
    result.data[1][3] = m->data[1][3];
    result.data[2][3] = m->data[2][3];
    result.data[3][3] = m->data[3][3];

    return result;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// lookat
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API fmat4 fmat4_lookat_agnostic(const float3 *eye, const float3 *target, const float3 *up)
{
    float3 sub = float3_sub(target, eye);
    float3 f = float3_normalize(&sub);
    float3 cross = float3_cross(&f, up);
    float3 s = float3_normalize(&cross);
    float3 u = float3_cross(&s, &f);
    fmat4 result = fmat4_identity();
    
    result.data[0][0] = s.xyz.x;
    result.data[0][1] = u.xyz.x;
    result.data[0][2] = -f.xyz.x;
    
    result.data[1][0] = s.xyz.y;
    result.data[1][1] = u.xyz.y;
    result.data[1][2] = -f.xyz.y;
    
    result.data[2][0] = s.xyz.z;
    result.data[2][1] = u.xyz.z;
    result.data[2][2] = -f.xyz.z;
    
    result.data[0][3] = -float3_dot(&s, eye);
    result.data[1][3] = -float3_dot(&u, eye);
    result.data[2][3] = float3_dot(&f, eye);
    
    return result;
}

VECMATH_API fmat4 fmat4_lookat_vulkan(const float3* eye, const float3* target, const float3* up)
{
    float3 sub = float3_sub(target, eye);
    float3 f = float3_normalize(&sub);
    float3 cross = float3_cross(up, &f);
    float3 r = float3_normalize(&cross);
    float3 u = float3_cross(&f, &r);
    
    fmat4 result = { 0 };
    
    result.matrix.m00 = r.xyz.x;
    result.matrix.m10 = r.xyz.y;
    result.matrix.m20 = r.xyz.z;
    result.matrix.m30 = -float3_dot(&r, eye);
    
    result.matrix.m01 = -u.xyz.x;  // Y flip for Vulkan
    result.matrix.m11 = -u.xyz.y;  // Y flip for Vulkan
    result.matrix.m21 = -u.xyz.z;  // Y flip for Vulkan
    result.matrix.m31 = float3_dot(&u, eye);
    
    result.matrix.m02 = -f.xyz.x;  // Right-handed: negative Z forward
    result.matrix.m12 = -f.xyz.y;  // Right-handed: negative Z forward
    result.matrix.m22 = -f.xyz.z;  // Right-handed: negative Z forward
    result.matrix.m32 = float3_dot(&f, eye);
    
    result.matrix.m03 = 0.0f;
    result.matrix.m13 = 0.0f;
    result.matrix.m23 = 0.0f;
    result.matrix.m33 = 1.0f;
    
    return result;
}

VECMATH_API fmat4 fmat4_lookat_directx(const float3* eye, const float3* target, const float3* up)
{
    float3 sub = float3_sub(target, eye);
    float3 f = float3_normalize(&sub);
    float3 cross0 = float3_cross(up, &f);
    float3 s = float3_normalize(&cross0);  // reverse cross for left-handed
    float3 cross1 = float3_cross(&f, &s);
    float3 u = float3_normalize(&cross1);  // reverse cross for left-handed
    
    fmat4 result = { 0 };
    result.data[0][0] = s.xyz.x;
    result.data[0][1] = u.xyz.x;
    result.data[0][2] = f.xyz.x;
    result.data[0][3] = 0.0f;
    
    result.data[1][0] = s.xyz.y;
    result.data[1][1] = u.xyz.y;
    result.data[1][2] = f.xyz.y;
    result.data[1][3] = 0.0f;
    
    result.data[2][0] = s.xyz.z;
    result.data[2][1] = u.xyz.z;
    result.data[2][2] = f.xyz.z;
    result.data[2][3] = 0.0f;
    
    result.data[3][0] = -float3_dot(&s, eye);
    result.data[3][1] = -float3_dot(&u, eye);
    result.data[3][2] = -float3_dot(&f, eye);
    result.data[3][3] = 1.0f;
    
    return result;
}

VECMATH_API fmat4 fmat4_lookat_opengl(const float3 *eye, const float3 *target, const float3 *up)
{
    float3 sub = float3_sub(target, eye);
    float3 f = float3_normalize(&sub);
    float3 cross0 = float3_cross(&f, up);
    float3 s = float3_normalize(&cross0);
    float3 cross1 = float3_cross(&s, &f);
    float3 u = float3_normalize(&cross1);
    
    fmat4 result = { 0 };
    
    result.data[0][0] = s.xyz.x;
    result.data[0][1] = u.xyz.x;
    result.data[0][2] = -f.xyz.x;  // negative for right-handed view space
    result.data[0][3] = 0.0f;
    
    result.data[1][0] = s.xyz.y;
    result.data[1][1] = u.xyz.y;
    result.data[1][2] = -f.xyz.y;  // negative for right-handed view space
    result.data[1][3] = 0.0f;
    
    result.data[2][0] = s.xyz.z;
    result.data[2][1] = u.xyz.z;
    result.data[2][2] = -f.xyz.z;  // negative for right-handed view space
    result.data[2][3] = 0.0f;
    
    result.data[3][0] = -float3_dot(&s, eye);
    result.data[3][1] = -float3_dot(&u, eye);
    result.data[3][2] = float3_dot(&f, eye);  // positive for right-handed
    result.data[3][3] = 1.0f;
    
    return result;
}

VECMATH_API dmat4 dmat4_lookat_agnostic(const double3 *eye, const double3 *target, const double3 *up)
{
    double3 sub = double3_sub(target, eye);
    double3 f = double3_normalize(&sub);
    double3 cross = double3_cross(&f, up);
    double3 s = double3_normalize(&cross);
    double3 u = double3_cross(&s, &f);
    dmat4 result = dmat4_identity();
    
    result.data[0][0] = s.xyz.x;
    result.data[0][1] = u.xyz.x;
    result.data[0][2] = -f.xyz.x;
    
    result.data[1][0] = s.xyz.y;
    result.data[1][1] = u.xyz.y;
    result.data[1][2] = -f.xyz.y;
    
    result.data[2][0] = s.xyz.z;
    result.data[2][1] = u.xyz.z;
    result.data[2][2] = -f.xyz.z;
    
    result.data[0][3] = -double3_dot(&s, eye);
    result.data[1][3] = -double3_dot(&u, eye);
    result.data[2][3] = double3_dot(&f, eye);
    
    return result;
}

VECMATH_API dmat4 dmat4_lookat_vulkan(const double3* eye, const double3* target, const double3* up)
{
    double3 sub = double3_sub(target, eye);
    double3 f = double3_normalize(&sub);
    double3 cross0 = double3_cross(&f, up);
    double3 s = double3_normalize(&cross0);
    double3 cross1 = double3_cross(&s, &f);
    double3 u = double3_normalize(&cross1);
    
    dmat4 result = { 0 };
    result.data[0][0] = s.xyz.x;
    result.data[0][1] = -s.xyz.y;  // flip Y for Vulkan
    result.data[0][2] = s.xyz.z;
    result.data[0][3] = -double3_dot(&s, eye);
    
    result.data[1][0] = u.xyz.x;
    result.data[1][1] = -u.xyz.y;  // flip Y for Vulkan
    result.data[1][2] = u.xyz.z;
    result.data[1][3] = -double3_dot(&u, eye);
    
    result.data[2][0] = -f.xyz.x;  // negative Z for right-handed view space
    result.data[2][1] = f.xyz.y;   // Y already flipped above
    result.data[2][2] = -f.xyz.z;  // negative Z for right-handed view space
    result.data[2][3] = double3_dot(&f, eye);  // positive for right-handed
    
    result.data[3][0] = 0.0f;
    result.data[3][1] = 0.0f;
    result.data[3][2] = 0.0f;
    result.data[3][3] = 1.0f;
    
    return result;
}

VECMATH_API dmat4 dmat4_lookat_directx(const double3* eye, const double3* target, const double3* up)
{
    double3 sub = double3_sub(target, eye);
    double3 f = double3_normalize(&sub);
    double3 cross0 = double3_cross(up, &f);
    double3 s = double3_normalize(&cross0);  // reverse cross for left-handed
    double3 cross1 = double3_cross(&f, &s);
    double3 u = double3_normalize(&cross1);  // reverse cross for left-handed
    
    dmat4 result = { 0 };
    result.data[0][0] = s.xyz.x;
    result.data[0][1] = u.xyz.x;
    result.data[0][2] = f.xyz.x;
    result.data[0][3] = 0.0f;
    
    result.data[1][0] = s.xyz.y;
    result.data[1][1] = u.xyz.y;
    result.data[1][2] = f.xyz.y;
    result.data[1][3] = 0.0f;
    
    result.data[2][0] = s.xyz.z;
    result.data[2][1] = u.xyz.z;
    result.data[2][2] = f.xyz.z;
    result.data[2][3] = 0.0f;
    
    result.data[3][0] = -double3_dot(&s, eye);
    result.data[3][1] = -double3_dot(&u, eye);
    result.data[3][2] = -double3_dot(&f, eye);
    result.data[3][3] = 1.0f;
    
    return result;
}

VECMATH_API dmat4 dmat4_lookat_opengl(const double3* eye, const double3* target, const double3* up)
{
    double3 sub = double3_sub(target, eye);
    double3 f = double3_normalize(&sub);
    double3 cross0 = double3_cross(&f, up);
    double3 s = double3_normalize(&cross0);
    double3 cross1 = double3_cross(&s, &f);
    double3 u = double3_normalize(&cross1);
    
    dmat4 result = { 0 };
    result.data[0][0] = s.xyz.x;
    result.data[0][1] = u.xyz.x;
    result.data[0][2] = -f.xyz.x;  // negative for right-handed view space
    result.data[0][3] = 0.0f;
    
    result.data[1][0] = s.xyz.y;
    result.data[1][1] = u.xyz.y;
    result.data[1][2] = -f.xyz.y;  // negative for right-handed view space
    result.data[1][3] = 0.0f;
    
    result.data[2][0] = s.xyz.z;
    result.data[2][1] = u.xyz.z;
    result.data[2][2] = -f.xyz.z;  // negative for right-handed view space
    result.data[2][3] = 0.0f;
    
    result.data[3][0] = -double3_dot(&s, eye);
    result.data[3][1] = -double3_dot(&u, eye);
    result.data[3][2] = double3_dot(&f, eye);  // positive for right-handed
    result.data[3][3] = 1.0f;
    
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// perspective
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API fmat4 fmat4_perspective_vulkan(float fov_rad, float aspect, float near, float far)
{
    float tan_half_fov = tanf(fov_rad * 0.5f);
    float f = 1.0f / tan_half_fov;
    float range_inv = 1.0f / (near - far);
    
    fmat4 result = { 0 };
    
    // column-major order for Vulkan
    result.data[0][0] = f / aspect;
    result.data[1][1] = f;
    result.data[2][2] = far * range_inv;
    result.data[2][3] = -1.0f;
    result.data[3][2] = far * near * range_inv;
    
    return result;
}

VECMATH_API fmat4 fmat4_perspective_directx(float fov_rad, float aspect, float near, float far)
{
    float tan_half_fov = tanf(fov_rad * 0.5f);
    float f = 1.0f / tan_half_fov;
    
    fmat4 result = { 0 };
    result.data[0][0] = f / aspect;
    result.data[1][1] = f;
    result.data[2][2] = far / (far - near);  // Z [0, 1] mapping
    result.data[2][3] = 1.0f;
    result.data[3][2] = -(far * near) / (far - near);
    
    return result;
}

VECMATH_API fmat4 fmat4_perspective_opengl(float fov_rad, float aspect, float near, float far)
{
    float tan_half_fov = tanf(fov_rad * 0.5f);
    float f = 1.0f / tan_half_fov;
    
    fmat4 result = { 0 };
    result.data[0][0] = f / aspect;
    result.data[1][1] = f;
    result.data[2][2] = (far + near) / (far - near);  // Z [-1, 1] mapping
    result.data[2][3] = 1.0f;
    result.data[3][2] = -(2.0f * far * near) / (far - near);
    
    return result;
}

VECMATH_API dmat4 dmat4_perspective_vulkan(double fov_rad, double aspect, double near, double far)
{
    double tan_half_fov = tan(fov_rad * 0.5f);
    double f = 1.0f / tan_half_fov;
    
    dmat4 result = { 0 };
    result.data[0][0] = f / aspect;
    result.data[1][1] = -f;  // flip Y for Vulkan's Y-down
    result.data[2][2] = far / (far - near);  // Z [0, 1] mapping
    result.data[2][3] = 1.0f;
    result.data[3][2] = -(far * near) / (far - near);
    
    return result;
}

VECMATH_API dmat4 dmat4_perspective_directx(double fov_rad, double aspect, double near, double far)
{
    double tan_half_fov = tan(fov_rad * 0.5f);
    double f = 1.0f / tan_half_fov;
    
    dmat4 result = { 0 };
    result.data[0][0] = f / aspect;
    result.data[1][1] = f;
    result.data[2][2] = far / (far - near);  // Z [0, 1] mapping
    result.data[2][3] = 1.0f;
    result.data[3][2] = -(far * near) / (far - near);
    
    return result;
}

VECMATH_API dmat4 dmat4_perspective_opengl(double fov_rad, double aspect, double near, double far)
{
    double tan_half_fov = tan(fov_rad * 0.5f);
    double f = 1.0f / tan_half_fov;
    
    dmat4 result = { 0 };
    result.data[0][0] = f / aspect;
    result.data[1][1] = f;
    result.data[2][2] = (far + near) / (far - near);  // Z [-1, 1] mapping
    result.data[2][3] = 1.0f;
    result.data[3][2] = -(2.0f * far * near) / (far - near);
    
    return result;
}

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// orthographic
/////////////////////////////////////////////////////////////////////////////////////

VECMATH_API fmat4 fmat4_orthographic_vulkan(float left, float right, float bottom, float top, float near, float far)
{
    float rl = right - left;
    float tb = top - bottom;
    float fn = far - near;
    
    fmat4 result = { 0 };
    result.data[0][0] = 2.0f / rl;
    result.data[1][1] = -2.0f / tb;  // flip Y for Vulkan
    result.data[2][2] = 1.0f / fn;   // Z [0, 1]
    result.data[3][0] = -(right + left) / rl;
    result.data[3][1] = (top + bottom) / tb;  // positive for Y-down
    result.data[3][2] = -near / fn;  // map near to 0, far to 1
    result.data[3][3] = 1.0f;
    
    return result;
}

VECMATH_API fmat4 fmat4_orthographic_directx(float left, float right, float bottom, float top, float near, float far)
{
    float rl = right - left;
    float tb = top - bottom;
    float fn = far - near;
    
    fmat4 result = { 0 };
    result.data[0][0] = 2.0f / rl;
    result.data[1][1] = 2.0f / tb;
    result.data[2][2] = 1.0f / fn;  // Z [0, 1]
    result.data[3][0] = -(right + left) / rl;
    result.data[3][1] = -(top + bottom) / tb;
    result.data[3][2] = -near / fn;  // map near to 0, far to 1
    result.data[3][3] = 1.0f;
    
    return result;
}

VECMATH_API fmat4 fmat4_orthographic_opengl(float left, float right, float bottom, float top, float near, float far)
{
    float rl = right - left;
    float tb = top - bottom;
    float fn = far - near;
    
    fmat4 result = { 0 };
    result.data[0][0] = 2.0f / rl;
    result.data[1][1] = 2.0f / tb;
    result.data[2][2] = -2.0f / fn;  // Z [-1, 1] (negative for right-handed)
    result.data[3][0] = -(right + left) / rl;
    result.data[3][1] = -(top + bottom) / tb;
    result.data[3][2] = -(far + near) / fn;  // map near to -1, far to 1
    result.data[3][3] = 1.0f;
    
    return result;
}

VECMATH_API dmat4 dmat4_orthographic_vulkan(double left, double right, double bottom, double top, double near, double far)
{
    double rl = right - left;
    double tb = top - bottom;
    double fn = far - near;
    
    dmat4 result = { 0 };
    result.data[0][0] = 2.0f / rl;
    result.data[1][1] = -2.0f / tb;  // flip Y for Vulkan
    result.data[2][2] = 1.0f / fn;   // Z [0, 1]
    result.data[3][0] = -(right + left) / rl;
    result.data[3][1] = (top + bottom) / tb;  // positive for Y-down
    result.data[3][2] = -near / fn;  // map near to 0, far to 1
    result.data[3][3] = 1.0f;
    
    return result;
}

VECMATH_API dmat4 dmat4_orthographic_directx(double left, double right, double bottom, double top, double near, double far)
{
    double rl = right - left;
    double tb = top - bottom;
    double fn = far - near;
    
    dmat4 result = { 0 };
    result.data[0][0] = 2.0f / rl;
    result.data[1][1] = 2.0f / tb;
    result.data[2][2] = 1.0f / fn;  // Z [0, 1]
    result.data[3][0] = -(right + left) / rl;
    result.data[3][1] = -(top + bottom) / tb;
    result.data[3][2] = -near / fn;  // map near to 0, far to 1
    result.data[3][3] = 1.0f;
    
    return result;
}

VECMATH_API dmat4 dmat4_orthographic_opengl(double left, double right, double bottom, double top, double near, double far)
{
    double rl = right - left;
    double tb = top - bottom;
    double fn = far - near;
    
    dmat4 result = { 0 };
    result.data[0][0] = 2.0f / rl;
    result.data[1][1] = 2.0f / tb;
    result.data[2][2] = -2.0f / fn;  // Z [-1, 1] (negative for right-handed)
    result.data[3][0] = -(right + left) / rl;
    result.data[3][1] = -(top + bottom) / tb;
    result.data[3][2] = -(far + near) / fn;  // map near to -1, far to 1
    result.data[3][3] = 1.0f;
    
    return result;
}

