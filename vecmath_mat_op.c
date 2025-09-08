#include "vecmath_mat_op.h"

/////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////// identity 
/////////////////////////////////////////////////////////////////////////////////////

fmat2 fmat2_identity()
{
    fmat2 mat;
    mat.matrix.m00 = 1.0f; mat.matrix.m01 = 0.0f;
    mat.matrix.m10 = 0.0f; mat.matrix.m11 = 1.0f;
    return mat;
}

fmat3 fmat3_identity()
{
    fmat3 mat;
    mat.matrix.m00 = 1.0f; mat.matrix.m01 = 0.0f; mat.matrix.m02 = 0.0f;
    mat.matrix.m10 = 0.0f; mat.matrix.m11 = 1.0f; mat.matrix.m12 = 0.0f;
    mat.matrix.m20 = 0.0f; mat.matrix.m21 = 0.0f; mat.matrix.m22 = 1.0f;
    return mat;
}

fmat4 fmat4_identity()
{
    fmat4 mat;
    mat.matrix.m00 = 1.0f; mat.matrix.m01 = 0.0f; mat.matrix.m02 = 0.0f; mat.matrix.m03 = 0.0f;
    mat.matrix.m10 = 0.0f; mat.matrix.m11 = 1.0f; mat.matrix.m12 = 0.0f; mat.matrix.m13 = 0.0f;
    mat.matrix.m20 = 0.0f; mat.matrix.m21 = 0.0f; mat.matrix.m22 = 1.0f; mat.matrix.m23 = 0.0f;
    mat.matrix.m30 = 0.0f; mat.matrix.m31 = 0.0f; mat.matrix.m32 = 0.0f; mat.matrix.m33 = 1.0f;
    return mat;
}

dmat2 dmat2_identity()
{
    dmat2 mat;
    mat.matrix.m00 = 1.0f; mat.matrix.m01 = 0.0f;
    mat.matrix.m10 = 0.0f; mat.matrix.m11 = 1.0f;
    return mat;
}

dmat3 dmat3_identity()
{
    dmat3 mat;
    mat.matrix.m00 = 1.0f; mat.matrix.m01 = 0.0f; mat.matrix.m02 = 0.0f;
    mat.matrix.m10 = 0.0f; mat.matrix.m11 = 1.0f; mat.matrix.m12 = 0.0f;
    mat.matrix.m20 = 0.0f; mat.matrix.m21 = 0.0f; mat.matrix.m22 = 1.0f;
    return mat;
}

dmat4 dmat4_identity()
{
    dmat4 mat;
    mat.matrix.m00 = 1.0f; mat.matrix.m01 = 0.0f; mat.matrix.m02 = 0.0f; mat.matrix.m03 = 0.0f;
    mat.matrix.m10 = 0.0f; mat.matrix.m11 = 1.0f; mat.matrix.m12 = 0.0f; mat.matrix.m13 = 0.0f;
    mat.matrix.m20 = 0.0f; mat.matrix.m21 = 0.0f; mat.matrix.m22 = 1.0f; mat.matrix.m23 = 0.0f;
    mat.matrix.m30 = 0.0f; mat.matrix.m31 = 0.0f; mat.matrix.m32 = 0.0f; mat.matrix.m33 = 1.0f;
    return mat;
}
