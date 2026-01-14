#ifndef VECMATH_MAT_OP_INCLUDED
#define VECMATH_MAT_OP_INCLUDED

#include "vecmath_defines.h"
#include "vecmath_types.h"

#ifdef __cplusplus 
extern "C" {
#endif

/// @brief returns the identity matrix
VECMATH_API fmat2 fmat2_identity();
VECMATH_API fmat3 fmat3_identity();
VECMATH_API fmat4 fmat4_identity();
VECMATH_API dmat2 dmat2_identity();
VECMATH_API dmat3 dmat3_identity();
VECMATH_API dmat4 dmat4_identity();

/// @brief multiplies a matrix with a vector
VECMATH_API fmat2 fmat2_mul_float2(const fmat2* m, const float2* v);
VECMATH_API fmat3 fmat3_mul_float3(const fmat3* m, const float3* v);
VECMATH_API fmat4 fmat4_mul_float4(const fmat4* m, const float4* v);
VECMATH_API dmat2 dmat2_mul_double2(const dmat2* m, const double2* v);
VECMATH_API dmat3 dmat3_mul_double3(const dmat3* m, const double3* v);
VECMATH_API dmat4 dmat4_mul_double4(const dmat4* m, const double4* v);

/// @brief transposes the matrix/ flip matrix by changing row and columns
VECMATH_API fmat2 fmat2_transpose(const fmat2* m);
VECMATH_API fmat3 fmat3_transpose(const fmat3* m); 
VECMATH_API fmat4 fmat4_transpose(const fmat4* m);
VECMATH_API dmat2 dmat2_transpose(const dmat2* m);
VECMATH_API dmat3 dmat3_transpose(const dmat3* m);
VECMATH_API dmat4 dmat4_transpose(const dmat4* m);

/// @brief calculates the determinant of a matrix
VECMATH_API float fmat2_determinant(const fmat2* m);
VECMATH_API float fmat3_determinant(const fmat3* m);
VECMATH_API float fmat4_determinant(const fmat4* m);
VECMATH_API double dmat2_determinant(const dmat2* m);
VECMATH_API double dmat3_determinant(const dmat3* m);
VECMATH_API double dmat4_determinant(const dmat4* m);

/// @brief calculates the matrix inverse
VECMATH_API fmat2 fmat2_inverse(const fmat2* m);
VECMATH_API fmat3 fmat3_inverse(const fmat3* m);
VECMATH_API fmat4 fmat4_inverse(const fmat4* m);
VECMATH_API dmat2 dmat2_inverse(const dmat2* m);
VECMATH_API dmat3 dmat3_inverse(const dmat3* m);
VECMATH_API dmat4 dmat4_inverse(const dmat4* m);

/// @brief matrix decomposition to retrieve translation
VECMATH_API float3 fmat4_get_translation_rowmajor(const fmat4* m);
VECMATH_API float3 fmat4_get_translation_colmajor(const fmat4* m);
VECMATH_API double3 dmat4_get_translation_rowmajor(const dmat4* m);
VECMATH_API double3 dmat4_get_translation_colmajor(const dmat4* m);

/// @brief matrix decomposition to retrieve scale
VECMATH_API float3 fmat4_get_scale_rowmajor(const fmat4* m);
VECMATH_API float3 fmat4_get_scale_colmajor(const fmat4* m);
VECMATH_API double3 dmat4_get_scale_rowmajor(const dmat4* m);
VECMATH_API double3 dmat4_get_scale_colmajor(const dmat4* m);

/// @brief matrix decomposition to retrieve rotation (uses euler angles)
VECMATH_API float3 fmat4_get_rotation_rowmajor(const fmat4* m);
VECMATH_API float3 fmat4_get_rotation_colmajor(const fmat4* m);
VECMATH_API double3 dmat4_get_rotation_rowmajor(const dmat4* m);
VECMATH_API double3 dmat4_get_rotation_colmajor(const dmat4* m);

/// @brief fully decomposes the matrix into it's components
VECMATH_API void fmat4_decompose_rowmajor(const fmat4* m, float3* translation, float3* rotation, float3* scale);
VECMATH_API void fmat4_decompose_colmajor(const fmat4* m, float3* translation, float3* rotation, float3* scale);
VECMATH_API void dmat4_decompose_rowmajor(const dmat4* m, double3* translation, double3* rotation, double3* scale);
VECMATH_API void dmat4_decompose_colmajor(const dmat4* m, double3* translation, double3* rotation, double3* scale);

/// @brief matrix translation
VECMATH_API fmat4 fmat4_translate_rowmajor(const fmat4* m, const float3* dir);
VECMATH_API fmat4 fmat4_translate_colmajor(const fmat4* m, const float3* dir);
VECMATH_API dmat4 dmat4_translate_rowmajor(const dmat4* m, const double3* dir);
VECMATH_API dmat4 dmat4_translate_colmajor(const dmat4* m, const double3* dir);

/// @brief matrix rotating
VECMATH_API fmat4 fmat4_rotate_colmajor(const fmat4* m, float angle, const float3* axis);
VECMATH_API fmat4 fmat4_rotate_rowmajor(const fmat4* m, float angle, const float3* axis);
VECMATH_API dmat4 dmat4_rotate_colmajor(const dmat4* m, double angle, const double3* axis);
VECMATH_API dmat4 dmat4_rotate_rowmajor(const dmat4 *m, double angle, const double3* axis);

/// @brief matrix scaleing
VECMATH_API fmat4 fmat4_scale_rowmajor(const fmat4* m, const float3* dim);
VECMATH_API fmat4 fmat4_scale_colmajor(const fmat4* m, const float3* dim);
VECMATH_API dmat4 dmat4_scale_rowmajor(const dmat4* m, const double3* dim);
VECMATH_API dmat4 dmat4_scale_colmajor(const dmat4* m, const double3* dim);

/// @brief lookat taking into account multiple rendering api coordinates
VECMATH_API fmat4 fmat4_lookat_agnostic(const float3* eye, const float3* target, const float3* up);
VECMATH_API fmat4 fmat4_lookat_vulkan(const float3* eye, const float3* target, const float3* up);
VECMATH_API fmat4 fmat4_lookat_directx(const float3* eye, const float3* target, const float3* up);
VECMATH_API fmat4 fmat4_lookat_opengl(const float3* eye, const float3* target, const float3* up);
VECMATH_API dmat4 dmat4_lookat_agnostic(const double3* eye, const double3* target, const double3* up);
VECMATH_API dmat4 dmat4_lookat_vulkan(const double3* eye, const double3* target, const double3* up);
VECMATH_API dmat4 dmat4_lookat_directx(const double3* eye, const double3* target, const double3* up);
VECMATH_API dmat4 dmat4_lookat_opengl(const double3* eye, const double3* target, const double3* up);

/// @brief perspective projection taking into account multiple rendering api coordinates
VECMATH_API fmat4 fmat4_perspective_vulkan(float fov_rad, float aspect, float nearVal, float farVal);
VECMATH_API fmat4 fmat4_perspective_directx(float fov_rad, float aspect, float nearVal, float farVal);
VECMATH_API fmat4 fmat4_perspective_opengl(float fov_rad, float aspect, float nearVal, float farVal);
VECMATH_API dmat4 dmat4_perspective_vulkan(double fov_rad, double aspect, double nearVal, double farVal);
VECMATH_API dmat4 dmat4_perspective_directx(double fov_rad, double aspect, double nearVal, double farVal);
VECMATH_API dmat4 dmat4_perspective_opengl(double fov_rad, double aspect, double nearVal, double farVal);

/// @brief orthographic projection taking into account multiple rendering api coordinates
VECMATH_API fmat4 fmat4_orthographic_vulkan(float left, float right, float bottom, float top, float near, float far);
VECMATH_API fmat4 fmat4_orthographic_directx(float left, float right, float bottom, float top, float near, float far);
VECMATH_API fmat4 fmat4_orthographic_opengl(float left, float right, float bottom, float top, float near, float far);
VECMATH_API dmat4 dmat4_orthographic_vulkan(double left, double right, double bottom, double top, double near, double far);
VECMATH_API dmat4 dmat4_orthographic_directx(double left, double right, double bottom, double top, double near, double far);
VECMATH_API dmat4 dmat4_orthographic_opengl(double left, double right, double bottom, double top, double near, double far);

#ifdef __cplusplus 
}
#endif

#endif // VECMATH_MAT_OP_INCLUDED