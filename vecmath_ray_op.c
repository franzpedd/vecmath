#include "vecmath_ray_op.h"
#include "vecmath_vec_op.h"
#include <math.h>

VECMATH_API fray fray_from_screen_point_vulkan(const float2* screenPos, const float2* windowSize, float fov, float aspectRatio, const float3* cameraPos, const float3* cameraFront, const float3* cameraUp)
{
    fray result;

    float ndcX = (2.0f * screenPos->xy.x) / windowSize->xy.x - 1.0f;
    float ndcY = 1.0f - (2.0f * screenPos->xy.y) / windowSize->xy.y;

    float3 cameraRight = float3_cross(cameraUp, cameraFront);
    cameraRight = float3_normalize(&cameraRight);

    float3 up = float3_cross(cameraFront, &cameraRight);
    float3 normalizedUp = float3_normalize(&up);

    float halfHeight = tanf(fov * 0.5f);
    float halfWidth = halfHeight * aspectRatio;

    float3 direction = {
        cameraFront->xyz.x + cameraRight.xyz.x * (ndcX * halfWidth) + normalizedUp.xyz.x * (ndcY * halfHeight),
        cameraFront->xyz.y + cameraRight.xyz.y * (ndcX * halfWidth) + normalizedUp.xyz.y * (ndcY * halfHeight),
        cameraFront->xyz.z + cameraRight.xyz.z * (ndcX * halfWidth) + normalizedUp.xyz.z * (ndcY * halfHeight)
    };

    result.origin = *cameraPos;
    result.direction = float3_normalize(&direction);

    return result;
}

VECMATH_API float3 fray_screen_to_world_point_vulkan(const float2* screenPos, const float2* windowSize, float distance, float fov, float aspectRatio, const float3* cameraPos, const float3* cameraFront, const float3* cameraUp)
{
    float ndcX = (2.0f * screenPos->xy.x) / windowSize->xy.x - 1.0f;
    float ndcY = 1.0f - (2.0f * screenPos->xy.y) / windowSize->xy.y;

    float3 aux = {
        cameraUp->xyz.y * cameraFront->xyz.z - cameraUp->xyz.z * cameraFront->xyz.y,
        cameraUp->xyz.z * cameraFront->xyz.x - cameraUp->xyz.x * cameraFront->xyz.z,
        cameraUp->xyz.x * cameraFront->xyz.y - cameraUp->xyz.y * cameraFront->xyz.x
    };

    float3 cameraRight = float3_normalize(&aux);

    float3 up = {
        cameraFront->xyz.y * cameraRight.xyz.z - cameraFront->xyz.z * cameraRight.xyz.y,
        cameraFront->xyz.z * cameraRight.xyz.x - cameraFront->xyz.x * cameraRight.xyz.z,
        cameraFront->xyz.x * cameraRight.xyz.y - cameraFront->xyz.y * cameraRight.xyz.x
    };
    float3 normalizedUp = float3_normalize(&up);

    float halfHeight = distance * tanf(fov * 0.5f);
    float halfWidth = halfHeight * aspectRatio;

    float3 worldPoint = {
        cameraPos->xyz.x + cameraFront->xyz.x * distance + cameraRight.xyz.x * (ndcX * halfWidth) + normalizedUp.xyz.x * (ndcY * halfHeight),
        cameraPos->xyz.y + cameraFront->xyz.y * distance + cameraRight.xyz.y * (ndcX * halfWidth) + normalizedUp.xyz.y * (ndcY * halfHeight),
        cameraPos->xyz.z + cameraFront->xyz.z * distance + cameraRight.xyz.z * (ndcX * halfWidth) + normalizedUp.xyz.z * (ndcY * halfHeight)
    };

    return worldPoint;
}