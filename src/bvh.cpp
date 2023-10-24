#include "bvh.h"
#include "draw.h"
#include "extra.h"
#include "interpolate.h"
#include "intersect.h"
#include "render.h"
#include "scene.h"
#include "texture.h"
#include <algorithm>
#include <bit>
#include <chrono>
#include <framework/opengl_includes.h>
#include <iostream>

#include <stack>
#include <cmath>
#include <glm/gtx/string_cast.hpp>
#include <queue>

int leafCreationHits = 0;
int dataNodeHits = 0;

// Helper method to fill in hitInfo object. This can be safely ignored (or extended).
// Note: many of the functions in this helper tie in to standard/extra features you will have
// to implement separately, see interpolate.h/.cpp for these parts of the project
void updateHitInfo(RenderState& state, const BVHInterface::Primitive& primitive, const Ray& ray, HitInfo& hitInfo)
{
    const auto& [v0, v1, v2] = std::tie(primitive.v0, primitive.v1, primitive.v2);
    const auto& mesh = state.scene.meshes[primitive.meshID];
    const auto n = glm::normalize(glm::cross(v1.position - v0.position, v2.position - v0.position));
    const auto p = ray.origin + ray.t * ray.direction;

    // First, fill in default data, unrelated to separate features
    hitInfo.material = mesh.material;
    hitInfo.normal = n;
    hitInfo.barycentricCoord = computeBarycentricCoord(v0.position, v1.position, v2.position, p);

    // Next, if `features.enableNormalMapping` is true, generate smoothly interpolated vertex normals
    if (state.features.enableNormalInterp) {
        hitInfo.normal = interpolateNormal(v0.normal, v1.normal, v2.normal, hitInfo.barycentricCoord);
    }

    // Next, if `features.enableTextureMapping` is true, generate smoothly interpolated vertex uvs
    if (state.features.enableTextureMapping) {
        hitInfo.texCoord = interpolateTexCoord(v0.texCoord, v1.texCoord, v2.texCoord, hitInfo.barycentricCoord);
    }

    // Finally, catch flipped normals
    if (glm::dot(ray.direction, n) > 0.0f) {
        hitInfo.normal = -hitInfo.normal;
    }
}

// BVH constructor; can be safely ignored. You should not have to touch this
// NOTE: this constructor is tested, so do not change the function signature.
BVH::BVH(const Scene& scene, const Features& features)
{
#ifndef NDEBUG
    // Store start of bvh build for timing
    using clock = std::chrono::high_resolution_clock;
    const auto start = clock::now();
#endif

    // Count the total nr. of triangles in the scene
    size_t numTriangles = 0;
    for (const auto& mesh : scene.meshes)
        numTriangles += mesh.triangles.size();

    // Given the input scene, gather all triangles over which to build the BVH as a list of Primitives
    std::vector<Primitive> primitives;
    primitives.reserve(numTriangles);
    for (uint32_t meshID = 0; meshID < scene.meshes.size(); meshID++) {
        const auto& mesh = scene.meshes[meshID];
        for (const auto& triangle : mesh.triangles) {
            primitives.push_back(Primitive {
                .meshID = meshID,
                .v0 = mesh.vertices[triangle.x],
                .v1 = mesh.vertices[triangle.y],
                .v2 = mesh.vertices[triangle.z] });
        }
    }

    // Tell underlying vectors how large they should approximately be
    m_primitives.reserve(numTriangles);
    m_nodes.reserve(numTriangles + 1);

    // Recursively build BVH structure; this is where your implementation comes in
    //std::cout << m_nodes.size();
    m_nodes.emplace_back(); // Create root node
    //std::cout << m_nodes.size();
    m_nodes.emplace_back(); // Create dummy node s.t. children are allocated on the same cache line
    buildRecursive(scene, features, primitives, RootIndex, 0);

    // Fill in boilerplate data
    buildNumLevels();
    buildNumLeaves();

#ifndef NDEBUG
    // Output end of bvh build for timing
    const auto end = clock::now();
    if (features.enableAccelStructure)
        std::cout << "acceleration enabled" << std::endl;
    std::cout << "BVH construction time: " << std::chrono::duration<double, std::milli>(end - start).count() << " ms, type: " << scene.type << ", lvl: " << m_numLevels << ", leaves: " << m_numLeaves << ", total triangles:" << m_primitives.size() << std::endl;
    std::cout << leafCreationHits << " triangles put in leaf nodes out of total" << std::endl;
    std::cout << dataNodeHits << " data nodes created" << std::endl;
#endif
}

// BVH helper method; allocates a new node and returns its index
// You should not have to touch this
uint32_t BVH::nextNodeIdx()
{
    const auto idx = static_cast<uint32_t>(m_nodes.size());
    m_nodes.emplace_back();
    return idx;
}

// TODO: Standard feature
// Given a BVH triangle, compute an axis-aligned bounding box around the primitive
// - primitive; a single triangle to be stored in the BVH
// - return;    an axis-aligned bounding box around the triangle
// This method is unit-tested, so do not change the function signature.
AxisAlignedBox computePrimitiveAABB(const BVHInterface::Primitive primitive)
{
    float x_min = glm::min(glm::min(primitive.v0.position[0], primitive.v1.position[0]), primitive.v2.position[0]);
    float y_min = glm::min(glm::min(primitive.v0.position[1], primitive.v1.position[1]), primitive.v2.position[1]);
    float z_min = glm::min(glm::min(primitive.v0.position[2], primitive.v1.position[2]), primitive.v2.position[2]);
    glm::vec3 lower { x_min, y_min, z_min };

    float x_max = glm::max(glm::max(primitive.v0.position[0], primitive.v1.position[0]), primitive.v2.position[0]);
    float y_max = glm::max(glm::max(primitive.v0.position[1], primitive.v1.position[1]), primitive.v2.position[1]);
    float z_max = glm::max(glm::max(primitive.v0.position[2], primitive.v1.position[2]), primitive.v2.position[2]);
    glm::vec3 upper { x_max, y_max, z_max };

   /* glm::vec3 lower = glm::min(glm::min(primitive.v0.position, primitive.v1.position), primitive.v2.position);
    glm::vec3 upper = glm::max(glm::max(primitive.v0.position, primitive.v1.position), primitive.v2.position);*/

    return { .lower = lower, .upper = upper };
}

// TODO: Standard feature
// Given a range of BVH triangles, compute an axis-aligned bounding box around the range.
// - primitive; a contiguous range of triangles to be stored in the BVH
// - return;    a single axis-aligned bounding box around the entire set of triangles
// This method is unit-tested, so do not change the function signature.
AxisAlignedBox computeSpanAABB(std::span<const BVHInterface::Primitive> primitives)
{
    float global_x_min = FLT_MAX, global_y_min = FLT_MAX, global_z_min = FLT_MAX;
    float global_x_max = -FLT_MAX, global_y_max = -FLT_MAX, global_z_max = -FLT_MAX;

    /*float global_x_min = primitives[0].v0.position[0], global_x_max = primitives[0].v0.position[0];
    float global_y_min = primitives[0].v0.position[1], global_y_max = primitives[0].v0.position[1];
    float global_z_min = primitives[0].v0.position[2], global_z_max = primitives[0].v0.position[2];*/
    /*if (primitives.size() == 1)
        return computePrimitiveAABB(primitives[0]);*/

    for (const auto& primitive : primitives) {
        float x_min = glm::min(glm::min(primitive.v0.position[0], primitive.v1.position[0]), primitive.v2.position[0]);
        float y_min = glm::min(glm::min(primitive.v0.position[1], primitive.v1.position[1]), primitive.v2.position[1]);
        float z_min = glm::min(glm::min(primitive.v0.position[2], primitive.v1.position[2]), primitive.v2.position[2]);

        float x_max = glm::max(glm::max(primitive.v0.position[0], primitive.v1.position[0]), primitive.v2.position[0]);
        float y_max = glm::max(glm::max(primitive.v0.position[1], primitive.v1.position[1]), primitive.v2.position[1]);
        float z_max = glm::max(glm::max(primitive.v0.position[2], primitive.v1.position[2]), primitive.v2.position[2]);

        global_x_min = glm::min(global_x_min, x_min);
        global_y_min = glm::min(global_y_min, y_min);
        global_z_min = glm::min(global_z_min, z_min);

        global_x_max = glm::max(global_x_max, x_max);
        global_y_max = glm::max(global_y_max, y_max);
        global_z_max = glm::max(global_z_max, z_max);
    }
    /*glm::vec3 lower = primitives[0].v0.position, upper = primitives[0].v0.position;
    for (const auto& primitive : primitives) {
        glm::vec3 lower = glm::min(glm::min(glm::min(primitive.v0.position, primitive.v1.position), primitive.v2.position),lower);
        glm::vec3 upper = glm::max(glm::max(glm::max(primitive.v0.position, primitive.v1.position), primitive.v2.position), upper);
    }*/

    /*for (const auto& primitive : primitives) {

    }*/

    glm::vec3 lower { global_x_min, global_y_min, global_z_min };
    glm::vec3 upper { global_x_max, global_y_max, global_z_max };

    return { .lower = lower, .upper = upper };
}

// TODO: Standard feature
// Given a BVH triangle, compute the geometric centroid of the triangle
// - primitive; a single triangle to be stored in the BVH
// - return;    the         geometric centroid of the triangle's vertices
// This method is unit-tested, so do not change the function signature.
glm::vec3 computePrimitiveCentroid(const BVHInterface::Primitive primitive)
{
    float x_center = (primitive.v0.position[0] + primitive.v1.position[0] + primitive.v2.position[0]) / 3.f;
    float y_center = (primitive.v0.position[1] + primitive.v1.position[1] + primitive.v2.position[1]) / 3.f;
    float z_center = (primitive.v0.position[2] + primitive.v1.position[2] + primitive.v2.position[2]) / 3.f;

    return glm::vec3 { x_center, y_center, z_center };
}

// TODO: Standard feature
// Given an axis-aligned bounding box, compute the longest axis; x = 0, y = 1, z = 2.
// - aabb;   the input axis-aligned bounding box
// - return; 0 for the x-axis, 1 for the y-axis, 2 for the z-axis
//           if several axes are equal in length, simply return the first of these
// This method is unit-tested, so do not change the function signature.
uint32_t computeAABBLongestAxis(const AxisAlignedBox& aabb)
{
    float x_distance = aabb.upper[0] - aabb.lower[0];
    float y_distance = aabb.upper[1] - aabb.lower[1];
    float z_distance = aabb.upper[2] - aabb.lower[2];

    if (x_distance >= y_distance && x_distance >= z_distance)
        return 0;
    if (y_distance >= z_distance)
        return 1;
    return 2;
}

int closerToLowerAABB(std::tuple<float, int> distance1, std::tuple<float, int> distance2)
{
    return std::get<0>(distance1) < std::get<0>(distance2);

    /*if (std::get<0>(distance1) < std::get<0>(distance2))
        return -1;
    if (std::get<0>(distance1) > std::get<0>(distance2))
        return 1;
    return 0;*/
}

void swapPrimitives(std::span<BVHInterface::Primitive> primitives, size_t index1, size_t index2)
{
    BVH::Primitive temp = primitives[index1];
    primitives[index1] = primitives[index2];
    primitives[index2] = temp;
}

int comparePrimitivesX(BVHInterface::Primitive p1, BVHInterface::Primitive p2) {
    return computePrimitiveCentroid(p1)[0] < computePrimitiveCentroid(p2)[0];
}
int comparePrimitivesY(BVHInterface::Primitive p1, BVHInterface::Primitive p2)
{
    return computePrimitiveCentroid(p1)[1] < computePrimitiveCentroid(p2)[1];
}
int comparePrimitivesZ(BVHInterface::Primitive p1, BVHInterface::Primitive p2)
{
    return computePrimitiveCentroid(p1)[2] < computePrimitiveCentroid(p2)[2];
}

// TODO: Standard feature
// Given a range of BVH triangles, sort these along a specified axis based on their geometric centroid.
// Then, find and return the split index in the range, such that the subrange containing the first element
// of the list is at least as big as the other, and both differ at most by one element in size.
// Hint: you should probably reuse `computePrimitiveCentroid()`
// - aabb;       the axis-aligned bounding box around the given triangle range
// - axis;       0, 1, or 2, determining on which axis (x, y, or z) the split must happen
// - primitives; the modifiable range of triangles that requires sorting/splitting along an axis
// - return;     the split position of the modified range of triangles
// This method is unit-tested, so do not change the function signature.
size_t splitPrimitivesByMedian(const AxisAlignedBox& aabb, uint32_t axis, std::span<BVHInterface::Primitive> primitives)
{
    using Primitive = BVHInterface::Primitive;
    int triangle_count = primitives.size();

    std::vector<std::tuple<float, int>> distances;
    distances.resize(triangle_count);

    // const auto& primitive : primitives
    for (int i = 0; i < triangle_count; i++) {
        const auto centroid = computePrimitiveCentroid(primitives[i]);
        /*std::cout << glm::to_string(primitives[i].v0.position) << ' '
                  << glm::to_string(primitives[i].v1.position) << ' '
                  << glm::to_string(primitives[i].v2.position) << " ===>>> "
                  << glm::to_string(centroid) << std::endl; */
        float distance = centroid[axis] - aabb.lower[axis];
        distances[i] = std::make_tuple(distance, i);
    }


    /*for (int i = 0; i < triangle_count; i++) 
        std::cout << std::get<0>(distances[i]) << " ";
    std::cout << " not updated 2 " << std::endl;*/

    //std::sort(distances.begin(), distances.end(), closerToLowerAABB);
    if (axis == 0)
        std::sort(primitives.begin(), primitives.end(), comparePrimitivesX);
    else if (axis == 1)
        std::sort(primitives.begin(), primitives.end(), comparePrimitivesY);
    else
        std::sort(primitives.begin(), primitives.end(), comparePrimitivesZ);
        
    /*for (int i = 0; i < triangle_count; i++)
        std::cout << std::get<0>(distances[i]) << " ";
     std::cout << " update 1 " << std::endl;*/
    //MIGHT BE WRONG +-1
    //this is fucked 
    //for (int i = 0; i < triangle_count / 2; i++) {
    //    //std::cout << std::get<0>(distances[i]) << " ";
    //    swapPrimitives(primitives, std::get<1>(distances[i]), i);
    //}
    /*std::cout << std::endl
              << std::endl;*/

    if (triangle_count % 2)
        return triangle_count / 2 + 1; // not sure if this returns the right index
    return triangle_count / 2;
}

// TODO: Standard feature
// Hierarchy traversal routine; called by the BVH's intersect(),
// you must implement this method and implement it carefully!
//
// If `features.enableAccelStructure` is not enabled, the method should just iterate the BVH's
// underlying primitives (or the scene's geometry). The default imlpementation already does this.
// You will have to implement the part which actually traverses the BVH for a faster intersect,
// given that `features.enableAccelStructure` is enabled.
//
// This method returns `true` if geometry was hit, and `false` otherwise. On first/closest hit, the
// distance `t` in the `ray` object is updated, and information is updated in the `hitInfo` object.
//
// - state;    the active scene, and a user-specified feature config object, encapsulated
// - bvh;      the actual bvh which should be traversed for faster intersection
// - ray;      the ray intersecting the scene's geometry
// - hitInfo;  the return object, with info regarding the hit geometry
// - return;   boolean, if geometry was hit or not
//
// This method is unit-tested, so do not change the function signature.
bool intersectRayWithBVH(RenderState& state, const BVHInterface& bvh, Ray& ray, HitInfo& hitInfo)
{
    // Relevant data in the constructed BVH
    std::span<const BVHInterface::Node> nodes = bvh.nodes();
    std::span<const BVHInterface::Primitive> primitives = bvh.primitives();

    using Node = BVHInterface::Node;
    using Primitive = BVHInterface::Primitive;
    // Return value
    bool is_hit = false;

    if (state.features.enableAccelStructure) {
    //if (false){
        // TODO: implement here your (probably stack-based) BVH traversal.
        //
        // Some hints (refer to bvh_interface.h either way). BVH nodes are packed, so the
        // data is not easily extracted. Helper methods are available, however:
        // - For a given node, you can test if the node is a leaf with `node.isLeaf()`.
        // - If the node is not a leaf, you can obtain the left/right children with `node.leftChild()` etc.
        // - If the node is a leaf, you can obtain the offset to and nr. of primitives in the bvh's list
        //   of underlying primitives with `node.primitiveOffset()` and `node.primitiveCount()`
        //
        // In short, you will have to step down the bvh, node by node, and intersect your ray
        // with the node's AABB. If this intersection passes, you should:
        // - if the node is a leaf, intersect with the leaf's primitives
        // - if the node is not a leaf, test the left and right children as well!
        //
        // Note that it is entirely possible for a ray to hit a leaf node, but not its primitives,
        // and it is likewise possible for a ray to hit both children of a node.
        /*for (Primitive primitive : primitives)
            std::cout << glm::to_string(primitive.v0.position) << "     " << glm::to_string(primitive.v1.position) << "     " << glm::to_string(primitive.v2.position) << "     " << std::endl;
        std::cout << "end here" << std::endl;*/
        std::stack<Node> stack;
        stack.push(nodes[0]);
        float t = ray.t;
        do {

            const Node curr = stack.top();
            stack.pop();
            
            
            //std::cout << "n: " << glm::to_string(hitInfo.normal) << "barCoord: " << glm::to_string(hitInfo.barycentricCoord);

            //std::cout << glm::to_string(curr.aabb.lower) << glm::to_string(curr.aabb.upper) << std::endl;
            // if we intersect with the node at the top
            if (intersectRayWithShape(curr.aabb, ray)) {
                ray.t = t;
               /* std::cout << "n: " << glm::to_string(hitInfo.normal) << "barCoord: " << glm::to_string(hitInfo.barycentricCoord);
                std::cout << std::endl
                          << std::endl;*/
               
            //if (true){
                //std::cout << curr.isLeaf() << " " << curr.primitiveOffset() << " " << curr.primitiveCount() << std::endl;
                // if we have a leaf node at the top
                if (curr.isLeaf()) {
                    //std::cout << "hit " << std::endl;
                    int offset = curr.primitiveOffset();
                    int count = curr.primitiveCount();
                    int final_i = offset + count;

                    for (int i = offset; i < final_i; i++) {
                        Primitive prim = primitives[i];
                        //Primitive prim = m_primitives[i];
                        const auto& [v0, v1, v2] = std::tie(prim.v0, prim.v1, prim.v2);
                        if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                            updateHitInfo(state, prim, ray, hitInfo);
                            is_hit = true;
                        }
                    }
                }

                // if we have a data node
                else {
                    //ray.t = std::numeric_limits<float>::max();
                    const Node left = nodes[curr.leftChild()];
                    const Node right = nodes[curr.rightChild()];

                    //if (intersectRayWithShape(left.aabb, ray))
                        stack.push(left);
                    //if (intersectRayWithShape(right.aabb, ray))
                        stack.push(right);
                }
            }

        } while (!stack.empty());

    } else {
        // Naive implementation; simply iterates over all primitives
        for (const auto& prim : primitives) {
            const auto& [v0, v1, v2] = std::tie(prim.v0, prim.v1, prim.v2);
            if (intersectRayWithTriangle(v0.position, v1.position, v2.position, ray, hitInfo)) {
                updateHitInfo(state, prim, ray, hitInfo);
                is_hit = true;
            }
        }
    }
    /*std::cout << std::endl
              << std::endl;*/
    // Intersect with spheres.
    for (const auto& sphere : state.scene.spheres)
        is_hit |= intersectRayWithShape(sphere, ray, hitInfo);

    return is_hit;
}

// TODO: Standard feature
// Leaf construction routine; you should reuse this in in `buildRecursive()`
// Given an axis-aligned bounding box, and a range of triangles, generate a valid leaf object
// and store the triangles in the `m_primitives` vector.
// You are free to modify this function's signature, as long as the constructor builds a BVH
// - scene;      the active scene
// - features;   the user-specified features object
// - aabb;       the axis-aligned bounding box around the primitives beneath this leaf
// - primitives; the range of triangles to be stored for this leaf
BVH::Node BVH::buildLeafData(const Scene& scene, const Features& features, const AxisAlignedBox& aabb, std::span<Primitive> primitives, uint32_t lastIndex)
{
    leafCreationHits += primitives.size();
    Node node;
    // TODO fill in the leaf's data; refer to `bvh_interface.h` for details

    if (features.enableAccelStructure) {
        // const auto& mesh_vertices = scene.meshes[primitives[0].meshID].vertices;

        /*auto it = std::find(mesh_vertices.begin(), mesh_vertices.end(), primitives[0]);
        int index = -1;
        if (it != mesh_vertices.end()) {
            index = std::distance(mesh_vertices.begin(), it);
        }*/

        //auto it = std::find(m_primitives.begin(), m_primitives.end(), primitives[0]);
        //int index = 0;
        ////if (it != m_primitives.end()) {
        //    index = std::distance(m_primitives.begin(), it);
        ////} else
        //   //std::cout << "2 possible poops above";

        int distance = std::distance(m_primitives.begin(), m_primitives.end());

        node.data[0] = 1u << 31; // leaf bit
        //this is wrong - lastIndex
        node.data[0] = node.data[0] | distance; // index value in the rest of the bits
        node.data[1] = uint32_t(primitives.size());
        node.aabb = computeSpanAABB(primitives);
    }

    // Copy the current set of primitives to the back of the primitives vector
    std::copy(primitives.begin(), primitives.end(), std::back_inserter(m_primitives));

    return node;
}

// TODO: Standard feature
// Node construction routine; you should reuse this in in `buildRecursive()`
// Given an axis-aligned bounding box, and left/right child indices, generate a valid node object.
// You are free to modify this function's signature, as long as the constructor builds a BVH
// - scene;           the active scene
// - features;        the user-specified features object
// - aabb;            the axis-aligned bounding box around the primitives beneath this node
// - leftChildIndex;  the index of the node's left child in `m_nodes`
// - rightChildIndex; the index of the node's right child in `m_nodes`
BVH::Node BVH::buildNodeData(const Scene& scene, const Features& features, const AxisAlignedBox& aabb, uint32_t leftChildIndex, uint32_t rightChildIndex)
{
    dataNodeHits++;
    Node node;
    // TODO fill in the node's data; refer to `bvh_interface.h` for details
    if (features.enableAccelStructure) {
        node.data[0] = uint32_t(leftChildIndex);
        node.data[1] = uint32_t(rightChildIndex);
        node.aabb = aabb;
    }
    return node;
}

// TODO: Standard feature
// Hierarchy construction routine; called by the BVH's constructor,
// you must implement this method and implement it carefully!
//
// You should implement the other BVH standard features first, and this feature last, as you can reuse
// most of the other methods to assemble this part. There are detailed instructions inside the
// method which we recommend you follow.
//
// Arguments:
// - scene;      the active scene
// - features;   the user-specified features object
// - primitives; a range of triangles to be stored in the BVH
// - nodeIndex;  index of the node you are currently working on, this is already allocated
//
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildRecursive(const Scene& scene, const Features& features, std::span<Primitive> primitives, uint32_t nodeIndex, uint32_t lastIndex)
{
    // WARNING: always use nodeIndex to index into the m_nodes array. never hold a reference/pointer,
    // because a push/emplace (in ANY recursive calls) might grow vectors, invalidating the pointers.

    // Compute the AABB of the current node.
    AxisAlignedBox aabb = computeSpanAABB(primitives);

    // As a starting point, we provide an implementation which creates a single leaf, and stores
    // all triangles inside it. You should remove or comment this, and work on your own recursive
    // construction algorithm that implements the following steps. Make sure to reuse the methods
    // you have previously implemented to simplify this process.
    //
    // 1. Determine if the node should be a leaf, when the nr. of triangles is less or equal to 4
    //    (hint; use the `LeafSize` constant)
    // 2. If it is a leaf, fill in the leaf's data, and store its range of triangles in `m_primitives`
    // 3. If it is a node:
    //    3a. Split the range of triangles along the longest axis into left and right subspans,
    //        using either median or SAH-Binning based on the `Features` object
    //    3b. Allocate left/right child nodes
    //        (hint: use `nextNodeIdx()`)
    //    3c. Fill in the current node's data; aabb, left/right child indices
    //    3d. Recursively build left/right child nodes over their respective triangles
    //        (hint; use `std::span::subspan()` to split into left/right ranges)

    // Just configure the current node as a giant leaf for now
    // m_nodes[nodeIndex] = buildLeafData(scene, features, aabb, primitives);
    if (features.enableAccelStructure) {
        if (primitives.size() <= BVH::LeafSize)
        //if (true)
            m_nodes[nodeIndex] = buildLeafData(scene, features, aabb, primitives, lastIndex);
        else {/*
            for (int i = 0; i < primitives.size(); i++)
                std::cout << glm::to_string(primitives[i].v0.position);*/
            int longAxis = computeAABBLongestAxis(aabb);
            //std::cout << longAxis;
            /*for (int i = 0; i < primitives.size(); i++)
                std::cout << computePrimitiveCentroid(primitives[i])[longAxis] - aabb.lower[longAxis] << " ";
            std::cout << "not updated  1" << std::endl;*/
            int split = splitPrimitivesByMedian(aabb, longAxis, primitives);
            //std::cout << split;
            /* for (int i = 0; i < primitives.size(); i++)
                std::cout << computePrimitiveCentroid(primitives[i])[longAxis] - aabb.lower[longAxis] << " ";
             std::cout << "updated  2 " << std::endl;*/
            int leftId = nextNodeIdx();
            int rightId = nextNodeIdx();

            auto leftSpan = primitives.subspan(0, split);
            auto rightSpan = primitives.subspan(split, std::dynamic_extent);//std::dynamic_extent
            m_nodes[nodeIndex] = buildNodeData(scene, features, aabb, leftId, rightId);
            buildRecursive(scene, features, leftSpan, leftId, lastIndex);
            buildRecursive(scene, features, rightSpan, rightId, lastIndex + leftSpan.size());
        }
    } else
        m_nodes[nodeIndex] = buildLeafData(scene, features, aabb, primitives, 0);
}

// TODO: Standard feature, or part of it
// Compute the nr. of levels in your hierarchy after construction; useful for `debugDrawLevel()`
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildNumLevels()
{
    m_numLevels = (uint32_t) floor(log2(m_nodes.size()));
    /*int count = 0;
    for (Node node : m_nodes)
        if (!node.isLeaf())
            count++;
    m_numLevels = count;*/
    
    /*using Node = BVHInterface::Node;
    std::queue<Node> nodes;
    nodes.push(m_nodes[0]);

    int dataCount = 0;
    do {
        Node curr = nodes.front();
        nodes.pop();
        if (! curr.isLeaf()){
            nodes.push(m_nodes[curr.leftChild()]);
            nodes.push(m_nodes[curr.rightChild()]);
            dataCount++;
        }
    } while (!nodes.empty());

    float levels = log2(dataCount + 1);
    assert( abs((int) levels - levels) > 0.001 );
    m_numLevels = (int)levels;*/
}

// Compute the nr. of leaves in your hierarchy after construction; useful for `debugDrawLeaf()`
// You are free to modify this function's signature, as long as the constructor builds a BVH
void BVH::buildNumLeaves()
{
    m_numLeaves = (uint32_t)(m_nodes.size() - pow(2, m_numLevels - 1) - 1);
    /*int count = 0;
    for (Node node : m_nodes)
        if (node.isLeaf())
            count++;
    m_numLeaves = count;*/
    /*using Node = BVHInterface::Node;
    std::queue<Node> nodes;
    nodes.push(m_nodes[0]);

    int leafCount = 0;
    do {
        Node curr = nodes.front();
        nodes.pop();
        if (curr.isLeaf())
            leafCount++;
        else {
            nodes.push(m_nodes[curr.leftChild()]);
            nodes.push(m_nodes[curr.rightChild()]);
        }
    } while (!nodes.empty());

    m_numLeaves = leafCount;*/
}
// Draw the bounding boxes of the nodes at the selected level. Use this function to visualize nodes
// for debugging. You may wish to implement `buildNumLevels()` first. We suggest drawing the AABB
// of all nodes on the selected level.
// You are free to modify this function's signature.
void BVH::debugDrawLevel(int level)
{
    // Example showing how to draw an AABB as a (white) wireframe box.
    // Hint: use draw functions (see `draw.h`) to draw the contained boxes with different
    // colors, transparencies, etc.
    //AxisAlignedBox aabb { .lower = glm::vec3(0.0f), .upper = glm::vec3(0.0f, 1.05f, 1.05f) };
    Sampler* sample = new Sampler(1);

    int firstOnLevel = pow(2, level) - 1;
    int lastOnLevel = pow(2, level + 1) - 1; // -2 in reality but use -1 bc for loop
    //if (level == m_numLevels)
    //    lastOnLevel = firstOnLevel + m_numLeaves;
    //for (int i = firstOnLevel; i < lastOnLevel; i++)
        //drawAABB(m_nodes[i].aabb, DrawMode::Wireframe, glm::vec3(sample->next_1d(), sample->next_1d(), sample->next_1d()), 1.f);
    //drawAABB(m_nodes[0].aabb, DrawMode::Wireframe, glm::vec3(sample->next_1d(), sample->next_1d(), sample->next_1d()), 1.f);
    
    /*for (int i = 0; i < m_nodes.size(); i++) {
        if (i == 1)
            continue;
        drawAABB(m_nodes[i].aabb, DrawMode::Wireframe, glm::vec3(sample->next_1d(), sample->next_1d(), sample->next_1d()), 1.f);
    }*/
    using Node = BVHInterface::Node;
    std::queue<Node> nodes;
    nodes.push(m_nodes[0]);

    if (level == 0) {
        drawAABB(m_nodes[0].aabb, DrawMode::Wireframe, glm::vec3(sample->next_1d(), sample->next_1d(), sample->next_1d()), 1.f);
        return;
    }
      
    int nodeIndex = 0;
    do {
        Node curr = nodes.front();
        nodes.pop();
        
        nodes.push(m_nodes[curr.leftChild()]);
        nodes.push(m_nodes[curr.rightChild()]);
    } while (!nodes.empty() && ++nodeIndex < firstOnLevel);

    while (!nodes.empty()) {
        Node curr = nodes.front();
        nodes.pop();

        drawAABB(curr.aabb, DrawMode::Wireframe, glm::vec3(sample->next_1d(), sample->next_1d(), sample->next_1d()), 1.f);
    }

}

// Draw data of the leaf at the selected index. Use this function to visualize leaf nodes
// for debugging. You may wish to implement `buildNumLeaves()` first. We suggest drawing the AABB
// of the selected leaf, and then its underlying primitives with different colors.
// - leafIndex; index of the selected leaf.
//              (Hint: not the index of the i-th node, but of the i-th leaf!)
// You are free to modify this function's signature.
void BVH::debugDrawLeaf(int leafIndex)
{
    // Example showing how to draw an AABB as a (white) wireframe box.
    // Hint: use drawTriangle (see `draw.h`) to draw the contained primitives
    //AxisAlignedBox aabb { .lower = glm::vec3(0.0f), .upper = glm::vec3(0.0f, 1.05f, 1.05f) };
    Sampler* sample = new Sampler(1);
    using Node = BVHInterface::Node;
    using Primitive = BVHInterface::Primitive;
    std::queue<Node> nodes;
    nodes.push(m_nodes[0]);

    int level = m_numLevels - 1;
    int firstOnLevel = pow(2, level) - 1;
    int lastOnLevel = pow(2, level + 1) - 1; // -2 in reality but use -1 bc for loop
    

    int nodeIndex = 0;
    if (!nodes.front().isLeaf()) {
        do {
            Node curr = nodes.front();
            nodes.pop();

            nodes.push(m_nodes[curr.leftChild()]);
            nodes.push(m_nodes[curr.rightChild()]);
        } while (!nodes.empty() && ++nodeIndex < firstOnLevel);
    }

    int counter = 0;
    while (!nodes.empty()) {
        Node curr = nodes.front();
        nodes.pop();

        if (counter <= leafIndex)
            for (int i = curr.primitiveOffset(); i < curr.primitiveOffset() + curr.primitiveCount(); i++) {
                Primitive p = m_primitives[i];
                drawTriangle(p.v0, p.v1, p.v2);
            }
        counter++;
    }
}