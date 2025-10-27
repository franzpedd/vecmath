#ifndef VECMATH_QUAT_OP_INCLUDED
#define VECMATH_QUAT_OP_INCLUDED

#include "vecmath_defines.h"
#include "vecmath_types.h"

#ifdef __cplusplus 
extern "C" {
#endif

/// @brief returns identity quaternion
VECMATH_API fquat fquat_identity();
VECMATH_API dquat dquat_identity();

/// @brief multiply two quaternions
VECMATH_API fquat fquat_mul(const fquat* q1, const fquat* q2);
VECMATH_API dquat dquat_mul(const dquat* q1, const dquat* q2);

/// @brief returns the quaternion length
VECMATH_API float fquat_length(const fquat* q);
VECMATH_API double dquat_length(const dquat* q);

/// @brief quaternion conjugate (inverse for unit quaternions)
VECMATH_API fquat fquat_conjugate(const fquat* q);
VECMATH_API dquat dquat_conjugate(const dquat* q);

/// @brief normalizes the quaternion
VECMATH_API fquat fquat_normalize(const fquat* q);
VECMATH_API dquat dquat_normalize(const dquat* q);

/// @brief returns the dot product of two quaternions
VECMATH_API float fquat_dot(const fquat* q1, const fquat* q2);
VECMATH_API double dquat_dot(const dquat* q1, const dquat* q2);

/// @brief performs linear interpolation
VECMATH_API fquat fquat_lerp(const fquat* q1, const fquat* q2, float t);
VECMATH_API dquat dquat_lerp(const dquat* q1, const dquat* q2, double t);

/// @brief performs spherical interpolation
VECMATH_API fquat fquat_slerp(const fquat* q1, const fquat* q2, float t);
VECMATH_API dquat dquat_slerp(const dquat* q1, const dquat* q2, double t);

// create quaternion from axis-angle representation
VECMATH_API fquat fquat_from_euler(const float3* rad);
VECMATH_API dquat dquat_from_euler(const double3* d);

// converts a quaternion into euler angles (row=x, pitch=y, yaw=z)
VECMATH_API float3 fquat_to_euler(const fquat* q);
VECMATH_API double3 dquat_to_euler(const dquat* q);

/// @brief converts a quaternion to a matrix
VECMATH_API fmat4 fquat_to_fmat4_rowmajor(const fquat* q);
VECMATH_API fmat4 fquat_to_fmat4_colmajor(const fquat* q);
VECMATH_API dmat4 dquat_to_dmat4_rowmajor(const dquat* q);
VECMATH_API dmat4 dquat_to_dmat4_colmajor(const dquat* q);

#ifdef __cplusplus 
}
#endif

#endif // VECMATH_QUAT_OP_INCLUDED