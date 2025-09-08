#ifndef VECMATH_MAT_OP_INCLUDED
#define VECMATH_MAT_OP_INCLUDED

#include "vecmath_types.h"

// identity, vec_mul, transpose, determinant, inverse, translation, rotation, scale, lookat, perspective, ortographic, decompose, get_translation, get_rotation, get_scale

/// @brief returns the identity matrix
fmat2 fmat2_identity();
fmat3 fmat3_identity();
fmat4 fmat4_identity();
dmat2 dmat2_identity();
dmat3 dmat3_identity();
dmat4 dmat4_identity();

/// @brief multiplies a matrix with a vector


#endif // VECMATH_MAT_OP_INCLUDED