#include "light.h"
#include "bvh_interface.h"
#include "config.h"
#include "draw.h"
#include "intersect.h"
#include "render.h"
#include "scene.h"
#include "shading.h"
#include "bvh.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
DISABLE_WARNINGS_POP()


// TODO: Standard feature
// Given a single segment light, transform a uniformly distributed 1d sample in [0, 1),
// into a uniformly sampled position and an interpolated color on the segment light,
// and write these into the reference return values.
// - sample;    a uniformly distributed 1d sample in [0, 1)
// - light;     the SegmentLight object, see `common.h`
// - position;  reference return value of the sampled position on the light
// - color;     reference return value of the color emitted by the light at the sampled position
// This method is unit-tested, so do not change the function signature.
void sampleSegmentLight(const float& sample, const SegmentLight& light, glm::vec3& position, glm::vec3& color)
{ 
    glm::vec3 start = light.endpoint0;
    glm::vec3 end = light.endpoint1;
    
    glm::vec3 startColor = light.color0; 
    glm::vec3 endColor = light.color1; 
    
    glm::vec3 direction = end - start;
    
    position = start + sample * direction;
    color = sample * endColor + (1 - sample) * startColor;
}

// TODO: Standard feature
// Given a single paralellogram light, transform a uniformly distributed 2d sample in [0, 1),
// into a uniformly sampled position and interpolated color on the paralellogram light,
// and write these into the reference return values.
// - sample;   a uniformly distributed 2d sample in [0, 1)
// - light;    the ParallelogramLight object, see `common.h`
// - position; reference return value of the sampled position on the light
// - color;    reference return value of the color emitted by the light at the sampled position
// This method is unit-tested, so do not change the function signature.
void sampleParallelogramLight(const glm::vec2& sample, const ParallelogramLight& light, glm::vec3& position, glm::vec3& color)
{
    glm::vec3 v0 = light.v0;
    glm::vec3 edge1 = light.edge01;
    glm::vec3 edge2 = light.edge02;

    float alpha = sample[0];
    float beta = sample[1];

    // New position is linear combination of the edge vectors
    position = v0 + alpha * edge1 + beta * edge2;

    // This is just the bilinear interpolation formula
    color = light.color0 * (1 - alpha) * (1 - beta)
            + light.color1 * alpha * (1 - beta)
            + light.color2 * (1 - alpha) * beta
            + light.color3 * alpha * beta;
}

// TODO: Standard feature
// Given a sampled position on some light, and the emitted color at this position, return whether
// or not the light is visible from the provided ray/intersection.
// For a description of the method's arguments, refer to 'light.cpp'
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        whether the light is visible (true) or not (false)
// This method is unit-tested, so do not change the function signature.
bool visibilityOfLightSampleBinary(RenderState& state, const glm::vec3& lightPosition, const glm::vec3 &lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    if (!state.features.enableShadows) {
        // Shadows are disabled in the renderer
        return true;
    } else {
        // Shadows are enabled in the renderer
        // TODO: implement this function; currently, the light simply passes through
        if (ray.t == std::numeric_limits<float>::max())
            return false;

        HitInfo temp = hitInfo;
        glm::vec3 intersectionPoint = ray.origin + ray.direction * ray.t;

        // Compute the direction from the intersection point to the light source.
        glm::vec3 toLight = lightPosition - intersectionPoint;
        float distanceToLight = glm::length(toLight);

        // Create a shadow ray from the intersection point to the light source with bias
        float bias = 0.0001f;

        Ray shadowRay;
        shadowRay.origin = intersectionPoint;
        shadowRay.direction = toLight / distanceToLight;

        // Offset the ray by the bias
        shadowRay.origin += shadowRay.direction * bias;

        // Check for intersections from shadowRay to lightsource
        intersectRayWithBVH(state, state.bvh, shadowRay, temp);
        if (shadowRay.t < distanceToLight) {
            return false;
        }

        // If no obstructions were found, the light is visible.
        return true;
    }
}

// TODO: Standard feature
// Given a sampled position on some light, and the emitted color at this position, return the actual
// light that is visible from the provided ray/intersection, or 0 if this is not the case.
// Use the following blending operation: lightColor = lightColor * kd * (1 - alpha)
// Please reflect within 50 words in your report on why this is incorrect, and illustrate
// two examples of what is incorrect.
//
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        the visible light color that reaches the intersection
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 visibilityOfLightSampleTransparency(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    glm::vec3 result = lightColor;
    glm::vec3 intersectionPoint = ray.origin + ray.t * ray.direction;

    glm::vec3 toLight = lightPosition - intersectionPoint;

    HitInfo temp = hitInfo;

    float bias = 0.0001f;
    Ray shadowRay;
    shadowRay.origin = intersectionPoint;
    shadowRay.direction = glm::normalize(toLight);
    shadowRay.origin += bias * shadowRay.direction;

    float distanceToLight = glm::length(toLight);

    // Iterate over the scene objects and check for intersections with the shadow ray.
    for (const auto& mesh : state.scene.meshes) {
        const std::vector<Vertex>& vertices = mesh.vertices;
        const std::vector<glm::uvec3>& triangles = mesh.triangles;

        for (const glm::uvec3& triangle : triangles) {
            const glm::vec3& vertex1 = vertices[triangle.x].position;
            const glm::vec3& vertex2 = vertices[triangle.y].position;
            const glm::vec3& vertex3 = vertices[triangle.z].position;

            // Perform ray-triangle intersection test.
            if (intersectRayWithTriangle(vertex1, vertex2, vertex3, shadowRay, temp)) {
                if (shadowRay.t < distanceToLight) {
                    float alpha = mesh.material.transparency;
                    result *= (1 - alpha);
                    break; // Stop checking other triangles for transparency.
                }
            }
        }
    }

    return result;
}

// TODO: Standard feature
// Given a single point light, compute its contribution towards an incident ray at an intersection point.
//
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the light is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;   the active scene, feature config, bvh, and a thread-safe sampler
// - light;   the PointLight object, see `common.h`
// - ray;     the incident ray to the current intersection
// - hitInfo; information about the current intersection
// - return;  reflected light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionPointLight(RenderState& state, const PointLight& light, const Ray& ray, const HitInfo& hitInfo)
{
    glm::vec3 color = visibilityOfLightSample(state, light.position, light.color, ray, hitInfo);


    // TODO: modify this function to incorporate visibility corerctly
    glm::vec3 p = ray.origin + ray.t * ray.direction;
    glm::vec3 l = glm::normalize(light.position - p);
    glm::vec3 v = -ray.direction;
    return computeShading(state, v, l, color, hitInfo);
}

// TODO: Standard feature
// Given a single segment light, compute its contribution towards an incident ray at an intersection point
// by integrating over the segment, taking `numSamples` samples from the light source.
//
// Hint: you can sample the light by using `sampleSegmentLight(state.sampler.next_1d(), ...);`, which
//       you should implement first.
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the sample is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;      the active scene, feature config, bvh, and a thread-safe sampler
// - light;      the SegmentLight object, see `common.h`
// - ray;        the incident ray to the current intersection
// - hitInfo;    information about the current intersection
// - numSamples; the number of samples you need to take
// - return;     accumulated light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionSegmentLight(RenderState& state, const SegmentLight& light, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples)
{
    // TODO: implement this function; repeat numSamples times:
    // - sample the segment light
    // - test the sample's visibility
    // - then evaluate the phong model
    glm::vec3 accumulatedLight(0.0f);
    glm::vec3 sampledPosition(0.0f);
    glm::vec3 sampledColor(0.0f);
    glm::vec3 color(0.0f);
    float N = numSamples;

    for (int i = 0; i < numSamples; i++) {

        // Sample the segment light 
        sampleSegmentLight(state.sampler.next_1d(), light, sampledPosition, sampledColor);

        // Compute the color
        color = visibilityOfLightSample(state, sampledPosition, sampledColor / N, ray, hitInfo);

        if (color == glm::vec3(0.0f))
            continue;

        glm::vec3 p = ray.origin + ray.t * ray.direction;
        glm::vec3 l = glm::normalize(sampledPosition - p);
        glm::vec3 v = -ray.direction;
        accumulatedLight += computeShading(state, v, l, color, hitInfo);
    }


    return accumulatedLight;
}

// TODO: Standard feature
// Given a single parralelogram light, compute its contribution towards an incident ray at an intersection point
// by integrating over the parralelogram, taking `numSamples` samples from the light source, and applying
// shading.
//
// Hint: you can sample the light by using `sampleParallelogramLight(state.sampler.next_1d(), ...);`, which
//       you should implement first.
// Hint: you should use `visibilityOfLightSample()` to account for shadows, and if the sample is visible, use
//       the result of `computeShading()`, whose submethods you should probably implement first in `shading.cpp`.
//
// - state;      the active scene, feature config, bvh, and a thread-safe sampler
// - light;      the ParallelogramLight object, see `common.h`
// - ray;        the incident ray to the current intersection
// - hitInfo;    information about the current intersection
// - numSamples; the number of samples you need to take
// - return;     accumulated light along the incident ray, based on `computeShading()`
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 computeContributionParallelogramLight(RenderState& state, const ParallelogramLight& light, const Ray& ray, const HitInfo& hitInfo, uint32_t numSamples)
{
    // TODO: implement this function; repeat numSamples times:
    // - sample the parallellogram light
    // - test the sample's visibility
    // - then evaluate the phong model
    glm::vec3 accumulatedLight(0.0f);
    glm::vec3 sampledPosition(0.0f);
    glm::vec3 sampledColor(0.0f);
    glm::vec3 color(0.0f);
    float N = numSamples;

    for (int i = 0; i < numSamples; i++) {

        // Sample the parallelogram light
        sampleParallelogramLight(state.sampler.next_2d(), light, sampledPosition, sampledColor);

        // Compute the color
        color = visibilityOfLightSample(state, sampledPosition, sampledColor / N, ray, hitInfo);

        glm::vec3 p = ray.origin + ray.t * ray.direction;
        glm::vec3 l = glm::normalize(sampledPosition - p);
        glm::vec3 v = -ray.direction;
        accumulatedLight += computeShading(state, v, l, color, hitInfo);
    }

    return accumulatedLight;
}

// This function is provided as-is. You do not have to implement it.
// Given a sampled position on some light, and the emitted color at this position, return the actual
// light that is visible from the provided ray/intersection, or 0 if this is not the case.
// This forowards to `visibilityOfLightSampleBinary`/`visibilityOfLightSampleTransparency` based on settings.
//
// - state;         the active scene, feature config, and the bvh
// - lightPosition; the sampled position on some light source
// - lightColor;    the sampled color emitted at lightPosition
// - ray;           the incident ray to the current intersection
// - hitInfo;       information about the current intersection
// - return;        the visible light color that reaches the intersection
//
// This method is unit-tested, so do not change the function signature.
glm::vec3 visibilityOfLightSample(RenderState& state, const glm::vec3& lightPosition, const glm::vec3& lightColor, const Ray& ray, const HitInfo& hitInfo)
{
    if (!state.features.enableShadows) {
        // Shadows are disabled in the renderer
        return lightColor;
    } else if (!state.features.enableTransparency) {
        // Shadows are enabled but transparency is disabled
        return visibilityOfLightSampleBinary(state, lightPosition, lightColor, ray, hitInfo) ? lightColor : glm::vec3(0);
    } else {
        // Shadows and transparency are enabled
        return visibilityOfLightSampleTransparency(state, lightPosition, lightColor, ray, hitInfo);
    }
}

// This function is provided as-is. You do not have to implement it.
glm::vec3 computeLightContribution(RenderState& state, const Ray& ray, const HitInfo& hitInfo)
{
    // Iterate over all lights
    glm::vec3 Lo { 0.0f };
    for (const auto& light : state.scene.lights) {
        if (std::holds_alternative<PointLight>(light)) {
            Lo += computeContributionPointLight(state, std::get<PointLight>(light), ray, hitInfo);
        } else if (std::holds_alternative<SegmentLight>(light)) {
            Lo += computeContributionSegmentLight(state, std::get<SegmentLight>(light), ray, hitInfo, state.features.numShadowSamples);
        } else if (std::holds_alternative<ParallelogramLight>(light)) {
            Lo += computeContributionParallelogramLight(state, std::get<ParallelogramLight>(light), ray, hitInfo, state.features.numShadowSamples);
        }
    }
    return Lo;
}