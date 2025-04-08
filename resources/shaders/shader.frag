#version 330 core
out vec4 fragColor;

// These are for the original FEM project
in vec4 normal_worldSpace;
in vec4 position_worldSpace;

uniform int wire = 0;
uniform float red = 1.0;
uniform float green = 1.0;
uniform float blue = 1.0;
uniform float alpha = 1.0;

// These are for our final project
in vec2 vUV;

uniform float gray; // For Voxel
uniform sampler3D densityTex; // For RayMarching
uniform mat4 invViewProj;

// Define the position of the box
const vec3 boxMin = vec3(-0.5);
const vec3 boxMax = vec3(0.5);

// Set the step parameters in ray marching
// If the simulation is not detailed enough, we shold decrease stepSize or increase maxSteps
uniform int size;
float voxelLength = 1.0 / float(size);

// Using AABB to detect the intersection
bool intersectBox(vec3 origin, vec3 dir, out float tMin, out float tMax) {
    vec3 invDir = 1.0 / dir;
    vec3 t0s = (boxMin - origin) * invDir;
    vec3 t1s = (boxMax - origin) * invDir;

    vec3 tSmalls = min(t0s, t1s);
    vec3 tBigs   = max(t0s, t1s);

    tMin = max(max(tSmalls.x, tSmalls.y), tSmalls.z);
    tMax = min(min(tBigs.x, tBigs.y), tBigs.z);

    return tMax >= tMin;
}

void main() {
    // Step1: vUV [0, 1] to [-1, 1] (NDC coord)
    vec2 screenPos = vUV * 2.0 - 1.0;

    // Step2: Construct near and far palnes' points
    vec4 pointNear = invViewProj * vec4(screenPos, 0.0, 1.0); // 0.0 for z-depth, and 1.0 just for homogenous coord, no special meaning here
    vec4 pointFar = invViewProj * vec4(screenPos, 1.0, 1.0); // Same above

    // Step3: Construct the ray
    vec3 rayOrigin = pointNear.xyz / pointNear.w; // Normalize it using our homogenous coord's value
    vec3 rayTarget = pointFar.xyz / pointFar.w;
    vec3 rayDir = normalize(rayTarget - rayOrigin);

    // the time that the ray enters the cube and exit the cube
    float tEnter, tExit;
    // if not intersect, discard this ray
    if (!intersectBox(rayOrigin, rayDir, tEnter, tExit)) discard;

    // Ray Marching Logic starts here
    vec4 finalColor = vec4(0.0);

    // What we want for steps: (tMax - tMin) / stepSize ≈ maxSteps //TODO:: Need a better strategy here
    float rayLength = tExit - tEnter;
    float stepSize = 0.5 / float(size); // half voxel
    int maxSteps = int(rayLength / stepSize) + 1;

    // Read the enter time and enter position
    float t = tEnter + 0.5 * stepSize;
    vec3 pos = rayOrigin + rayDir * t;

    for (int i = 0; i < maxSteps && t < tExit; ++i) {
        // Map the pos to [0, 1] to get the related density because density is converted into [0, 1]^3 3D image
        vec3 texCoord = (pos - boxMin) / (boxMax - boxMin);

        // if texCoord out of boundaries, exit
        if (any(lessThan(texCoord, vec3(0.0))) || any(greaterThan(texCoord, vec3(1.0))))
            break;

        // We can then get the density
        float density = texture(densityTex, texCoord).r;

        float alpha = density;
        vec3 color = vec3(1-density);

        // Alpha Blending
        finalColor.rgb += (1.0 - finalColor.a) * alpha * color;
        finalColor.a += (1.0 - finalColor.a) * alpha;

        // Exit earlier if the alpha is already large enough
        if (finalColor.a > 0.95)
            break;

        t += stepSize;
        pos += rayDir * stepSize;
    }

    if (finalColor.a == 0.0)
        discard;

    fragColor = clamp(finalColor, 0.0, 1.0);
}
