#ifndef VECMATH_RAY_OP_INCLUDED
#define VECMATH_RAY_OP_INCLUDED

#include "vecmath_defines.h"
#include "vecmath_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/// @brief  casts a ray from screen coordinate 2d point (usually mouse coordinates)
VECMATH_API fray fray_from_screen_point_vulkan(const float2* screenPos, const float2* viewportSize, const fmat4* inverseProj, const fmat4* inverseView);

/// @brief returns a 3d world point based on ray and distance
VECMATH_API float3 fray_get_point(const fray* ray, float distance);

/// @brief checks if the ray intersects with different objects
VECMATH_API vecbool fray_sphere_intersection(const fray* ray, const float3* center, float radius, float* outDistance);
VECMATH_API vecbool fray_quad_intersection(const fray* ray, const float3* v0, const float3* v1, const float3* v2, const float3* v3, float* outDistance);
VECMATH_API vecbool fray_aabb_intersection(const fray* ray, const float3* minBound, const float3* maxBound, float* outDistance);
VECMATH_API vecbool fray_obb_intersection(const fray* ray, const float3* center, const fmat3* rotation, const float3* halfExtents, float* outDistance);

#ifdef __cplusplus
}
#endif

#endif // VECMATH_RAY_INCLUDED