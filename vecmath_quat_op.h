#ifndef VECMATH_QUAT_OP_INCLUDED
#define VECMATH_QUAT_OP_INCLUDED

#include "vecmath_types.h"

/// @brief returns identity quaternion
fquat fquat_identity();
dquat dquat_identity();

/// @brief multiply two quaternions
fquat fquat_mul(const fquat* q1, const fquat* q2);
dquat dquat_mul(const dquat* q1, const dquat* q2);

/// @brief returns the quaternion length
float fquat_length(const fquat* q);
double dquat_length(const dquat* q);

/// @brief quaternion conjugate (inverse for unit quaternions)
fquat fquat_conjugate(const fquat* q);
dquat dquat_conjugate(const dquat* q);

/// @brief normalizes the quaternion
fquat fquat_normalize(const fquat* q);
dquat dquat_normalize(const dquat* q);

/// @brief returns the dot product of two quaternions
float fquat_dot(const fquat* q1, const fquat* q2);
double dquat_dot(const dquat* q1, const dquat* q2);

/// @brief performs linear interpolation
fquat fquat_lerp(const fquat* q1, const fquat* q2, float t);
dquat dquat_lerp(const dquat* q1, const dquat* q2, double t);

/// @brief performs spherical interpolation
fquat fquat_slerp(const fquat* q1, const fquat* q2, float t);
dquat dquat_slerp(const dquat* q1, const dquat* q2, double t);

// create quaternion from axis-angle representation
fquat fquat_from_euler(const float3* axis, float angle_rad);
dquat dquat_from_euler(const double3* axis, double angle_rad);

// converts a quaternion into euler angles (row=x, pitch=y, yaw=z)
float3 fquat_to_euler(const fquat* q);
double3 dquat_to_euler(const dquat* q);

/// @brief converts a quaternion to a matrix
fmat4 fquat_to_fmat4_rowmajor(const fquat* q);
fmat4 fquat_to_fmat4_colmajor(const fquat* q);
dmat4 dquat_to_dmat4_rowmajor(const dquat* q);
dmat4 dquat_to_dmat4_colmajor(const dquat* q);

#endif // VECMATH_QUAT_OP_INCLUDED