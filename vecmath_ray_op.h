#ifndef VECMATH_RAY_OP_INCLUDED
#define VECMATH_RAY_OP_INCLUDED

#include "vecmath_defines.h"
#include "vecmath_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/// @brief casts a ray from screen coordinate 2d point (usually mouse coordinates)
VECMATH_API fray fray_from_screen_point_vulkan(const float2* screenPos, const float2* windowSize, float fov, float aspectRatio, const float3* cameraPos, const float3* cameraFront, const float3* cameraUp);

/// @brief returns the 3D point at a distance from a screen position (usually mouse coordinates)
VECMATH_API float3 fray_screen_to_world_point_vulkan(const float2* screenPos, const float2* windowSize, float distance, float fov, float aspectRatio, const float3* cameraPos, const float3* cameraFront, const float3* cameraUp);

#ifdef __cplusplus
}
#endif

#endif // VECMATH_RAY_INCLUDED