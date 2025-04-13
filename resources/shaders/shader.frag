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

// Fluid style: 0=gray, 1=blue water
uniform int fluidStyle = 1;  // Default to blue water

// Define the position of the box
const vec3 boxMin = vec3(-0.5);
const vec3 boxMax = vec3(0.5);

// Set the step parameters in ray marching
// If the simulation is not detailed enough, we should decrease stepSize or increase maxSteps
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

// Helper functions for water effects
float fresnel(float cosTheta, float F0) {
    return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}

void main() {
    // Step1: vUV [0, 1] to [-1, 1] (NDC coord)
    vec2 screenPos = vUV * 2.0 - 1.0;

    // Step2: Construct near and far planes' points
    vec4 pointNear = invViewProj * vec4(screenPos, 0.0, 1.0);
    vec4 pointFar = invViewProj * vec4(screenPos, 1.0, 1.0);

    // Step3: Construct the ray
    vec3 rayOrigin = pointNear.xyz / pointNear.w;
    vec3 rayTarget = pointFar.xyz / pointFar.w;
    vec3 rayDir = normalize(rayTarget - rayOrigin);

    // The time that the ray enters the cube and exits the cube
    float tEnter, tExit;
    // If not intersect, discard this ray
    if (!intersectBox(rayOrigin, rayDir, tEnter, tExit)) discard;

    // Ray Marching Logic starts here
    vec4 finalColor = vec4(0.0);

    // What we want for steps: (tMax - tMin) / stepSize ≈ maxSteps
    float rayLength = tExit - tEnter;
    float stepSize = 0.5 / float(size);  // Half voxel
    int maxSteps = int(rayLength / stepSize) + 1;

    // Read the enter time and enter position
    float t = tEnter + 0.5 * stepSize;
    vec3 pos = rayOrigin + rayDir * t;
    
    // For edge highlighting
    float prevDensity = 0.0;
    vec3 prevPos = pos;

    // Add overall water colors
    vec3 waterBaseColor = vec3(0.1, 0.4, 0.8);  // Base blue
    vec3 waterDeepColor = vec3(0.0, 0.1, 0.3);  // Deep water blue
    vec3 waterHighlightColor = vec3(0.8, 0.9, 1.0);  // White foam/highlight
    vec3 waterSurfaceColor = vec3(0.5, 0.7, 1.0);  // Surface blue
    
    
    // For water refractions
    vec3 lightDir = normalize(vec3(0.5, 1.0, 0.5));
    float NdotL = max(dot(rayDir, lightDir), 0.0);

    for (int i = 0; i < maxSteps && t < tExit; ++i) {
        // Map the pos to [0, 1] to get the related density
        vec3 texCoord = (pos - boxMin) / (boxMax - boxMin);

        // If texCoord out of boundaries, exit
        if (any(lessThan(texCoord, vec3(0.0))) || any(greaterThan(texCoord, vec3(1.0))))
            break;

        // Get the density
        float density = texture(densityTex, texCoord).r;

        // Skip low density regions for efficiency
        if (density < 0.01) {
            t += stepSize * 2.0;
            pos += rayDir * stepSize * 2.0;
            continue;
        }
        
        // Calculate edge factor (density gradient)
        float edge = 0.0;
        if (i > 0) {
            edge = abs(density - prevDensity) * 15.0;
        }
        prevDensity = density;
        prevPos = pos;

        // Set color and alpha based on fluid style
        vec3 color;
        float alphaValue;

        if (fluidStyle == 0) {  // Gray (original)
            color = vec3(1.0 - density);
            alphaValue = density * 0.15;
        }
        else if (fluidStyle == 1) {  // Blue water
            // Base water color based on depth
            float depth = texCoord.y;
            color = mix(waterDeepColor, waterBaseColor, depth);
            
            // Add foam effect to high density areas
            color = mix(color, waterHighlightColor, pow(density, 3.0) * 0.5);
            
            // Add surface reflection at top of fluid, smooth transition that affects the top 20% of the fluid
            float surfaceFactor = smoothstep(0.8, 1.0, texCoord.y);
            color = mix(color, waterSurfaceColor, surfaceFactor * 0.5); //light blue
            
            // Add edge highlights for whose density changes rapidly
            color = mix(color, waterHighlightColor, clamp(edge, 0.0, 0.5));
            
            alphaValue = density * 0.15 * (1.0 + 0.5 * depth);
            
            // Add Fresnel effect at glancing angles (path)
            float fresnelEffect = fresnel(abs(dot(rayDir, vec3(0.0, 1.0, 0.0))), 0.02);
            color = mix(color, waterSurfaceColor, fresnelEffect * 0.3);
        }
       

        // Alpha Blending
        finalColor.rgb += (1.0 - finalColor.a) * alphaValue * color;
        finalColor.a += (1.0 - finalColor.a) * alphaValue;

        // Exit earlier if the alpha is already large enough
        if (finalColor.a > 0.95)
            break;

        t += stepSize;
        pos += rayDir * stepSize;
    }

    // Post-processing effects based on fluid style
    if (finalColor.a > 0.0) {
        if (fluidStyle == 1) {  // Blue water
            // Add subtle blue tint to the entire scene
            finalColor.rgb = mix(finalColor.rgb, vec3(0.4, 0.6, 0.9), 0.05);
            
            // Add subtle caustics effect
            float caustic = sin(vUV.x * 20.0) * sin(vUV.y * 20.0) * 0.5 + 0.5;
            finalColor.rgb += vec3(0.0, 0.0, 0.2) * caustic * 0.05;
        }
       
        
        // Add subtle vignette effect for all styles
        float vignette = 1.0 - dot(vUV - 0.5, vUV - 0.5) * 0.5;
        finalColor.rgb *= vignette;
        
        // Enhance contrast slightly
        finalColor.rgb = pow(finalColor.rgb, vec3(1.1));
    }

    if (finalColor.a == 0.0)
        discard;

    fragColor = clamp(finalColor, 0.0, 1.0);
}