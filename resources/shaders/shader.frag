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
uniform sampler3D colorTex;   // For per-cell color mapping
uniform mat4 invViewProj;
uniform float time; // For time-based animation effects

uniform int fluidStyle;

// Color map types - extended with new gradient effects
uniform int colorMapType = 0; // 0: default, 1: blue, 2: purple, 3: cyan-yellow, 4: orange-grey
                             // 5: fire gradient, 6: ocean gradient, 7: plasma gradient
                             // 8: rainbow, 9: aurora, 10: lava

uniform int renderMode = 0; // 0: volume, 1: shell

// Define the position of the box - increased from 0.5 to 2.0 for larger display
const vec3 boxMin = vec3(-2.0);
const vec3 boxMax = vec3(2.0);

// Set the step parameters in ray marching
uniform int size;
float voxelLength = 1.0 / float(size);

// Noise functions for visual effects
float hash(float n) {
    return fract(sin(n) * 43758.5453);
}

float noise(vec3 x) {
    vec3 p = floor(x);
    vec3 f = fract(x);
    f = f * f * (3.0 - 2.0 * f);

    float n = p.x + p.y * 157.0 + 113.0 * p.z;
    return mix(mix(mix(hash(n + 0.0), hash(n + 1.0), f.x),
                   mix(hash(n + 157.0), hash(n + 158.0), f.x), f.y),
               mix(mix(hash(n + 113.0), hash(n + 114.0), f.x),
                   mix(hash(n + 270.0), hash(n + 271.0), f.x), f.y), f.z);
}

float fbm(vec3 x) {
    float v = 0.0;
    float a = 0.5;
    vec3 shift = vec3(100);

    // Rotate to reduce axial bias
    mat3 rot = mat3(cos(0.5), sin(0.5), 0,
                    -sin(0.5), cos(0.5), 0,
                    0, 0, 1);

    for (int i = 0; i < 4; ++i) {
        v += a * noise(x);
        x = rot * x * 2.0 + shift;
        a *= 0.5;
    }
    return v;
}

// Fresnel effect calculation
float fresnel(float cosTheta, float F0) {
    return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
}

// Ray-box intersection test
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

// Get fluid color based on type, depth, and other parameters
vec3 getFluidColor(int effectiveColorType, float depth, float density, float edge) {
    vec3 color;
    vec3 waterBaseColor = vec3(0.1, 0.4, 0.8);  // Base blue
    vec3 waterDeepColor = vec3(0.0, 0.1, 0.3);  // Deep water blue
    vec3 waterHighlightColor = vec3(0.8, 0.9, 1.0);  // White foam/highlight
    vec3 waterSurfaceColor = vec3(0.5, 0.7, 1.0);  // Surface blue

    // Default water color (type 0)
    if (effectiveColorType == 0) {
        // Original blue water color logic with enhanced depth effect
        float depthEffect = pow(depth, 1.2); // Makes depth effect more pronounced
        color = mix(waterDeepColor, waterBaseColor, depthEffect);

        // Add foam effect to high density areas
        color = mix(color, waterHighlightColor, pow(density, 2.5) * 0.6);

        // Add surface reflection at top of fluid with stronger effect
        float surfaceFactor = smoothstep(0.75, 1.0, depth);
        color = mix(color, waterSurfaceColor, surfaceFactor * 0.6);

        // Add edge highlights with more intensity
        color = mix(color, waterHighlightColor, clamp(edge * 1.2, 0.0, 0.6));
    }
    // Cyan-blue (type 1)
    else if (effectiveColorType == 1) {
        // Cyan-blue theme with enhanced contrast
        color = mix(vec3(0.0, 0.9, 1.0), vec3(0.0, 0.1, 0.6), pow(1.0-depth, 1.5));
        // Add highlights at density variations
        color = mix(color, vec3(0.6, 0.9, 1.0), pow(density, 3.0) * 0.7);
    }
    // Purple electric (type 2)
    else if (effectiveColorType == 2) {
        // Purple electric theme with more dramatic gradient
        color = mix(vec3(0.7, 0.0, 1.0), vec3(0.1, 0.0, 0.4), pow(1.0-depth, 1.3));
        // Add highlights
        color = mix(color, vec3(0.9, 0.5, 1.0), pow(density, 2.8) * 0.6);
    }
    // Cyan-yellow (type 3)
    else if (effectiveColorType == 3) {
        // Cyan-yellow gradient with stronger contrast
        float gradientPos = depth * 0.7 + density * 0.3;
        color = mix(vec3(0.0, 0.9, 0.9), vec3(1.0, 0.9, 0.1), pow(gradientPos, 1.2));
    }
    // Orange-grey (type 4)
    else if (effectiveColorType == 4) {
        // Orange-grey explosion effect with more contrast
        float gradientPos = depth * 0.6 + density * 0.4;
        color = mix(vec3(1.0, 0.4, 0.0), vec3(0.6, 0.6, 0.6), pow(gradientPos, 1.1));
    }
    // Fire (type 5)
    else if (effectiveColorType == 5) {
        // Create a noise-based distortion of the depth
        float noiseScale = 5.0;
        float noiseStrength = 0.2;

        // Generate base noise - varies in all three dimensions
        float baseNoise = fbm(vec3(depth * noiseScale,
                                   density * noiseScale * 0.5,
                                   time * 0.2));

        // Create more detailed noise for fine details
        float detailNoise = noise(vec3(depth * noiseScale * 2.0,
                                      density * noiseScale * 3.0,
                                      time * 0.5));

        // Combine noise with depth to create irregular layering
        float layerPos = depth + (baseNoise * 0.3 + detailNoise * 0.1) * noiseStrength;
        // Add density influence - denser regions burn hotter (more yellow/white)
        layerPos += density * 0.15;

        // Create natural fire color gradient
        if (layerPos < 0.25) {
            // Inner core - bright white-yellow
            color = mix(vec3(1.0, 0.9, 0.5), vec3(1.0, 0.8, 0.0), smoothstep(0.0, 0.25, layerPos));
        } else if (layerPos < 0.5) {
            // Mid flame - yellow to orange
            color = mix(vec3(1.0, 0.8, 0.0), vec3(1.0, 0.5, 0.0), smoothstep(0.25, 0.5, layerPos));
        } else if (layerPos < 0.75) {
            // Outer flame - orange to red
            color = mix(vec3(1.0, 0.5, 0.0), vec3(0.9, 0.1, 0.0), smoothstep(0.5, 0.75, layerPos));
        } else {
            // Fading edges - red to dark red
            color = mix(vec3(0.9, 0.1, 0.0), vec3(0.4, 0.0, 0.0), smoothstep(0.75, 1.0, layerPos));
        }

        // Add flickering effect
        float flicker = noise(vec3(depth, density, time * 2.0)) * 0.2 + 0.9;
        color *= flicker;

        // Add glow based on density and position
        float glowStrength = 0.4;
        vec3 glowColor = vec3(1.0, 0.3, 0.0) * density * glowStrength;
        color += glowColor * (1.0 - layerPos);
    }
    // Ocean (type 6)
    else if (effectiveColorType == 6) {
        // Create noise-based distortion of the depth with increased strength
        float noiseScale = 4.0;
        float noiseStrength = 0.35;

        // Generate base noise for large flowing patterns
        float baseNoise = fbm(vec3(depth * noiseScale,
                                  density * noiseScale * 0.3,
                                  time * 0.3));

        // Create more detailed noise for small currents and eddies
        float detailNoise = noise(vec3(depth * noiseScale * 3.0,
                                      density * noiseScale * 2.0,
                                      time * 0.5));

        // Combine noise with depth for irregular layering
        float layerPos = depth + (baseNoise * 0.5 + detailNoise * 0.2) * noiseStrength;
        // Add slight density influence - denser regions appear deeper
        layerPos -= density * 0.15;

        // Create more vibrant ocean color gradient
        if (layerPos < 0.25) {
            // Deep ocean - much darker blue to create stronger contrast
            color = vec3(0.0, 0.0, 0.2);
        } else if (layerPos < 0.5) {
            // Mid-depths - rich royal blue
            color = mix(vec3(0.0, 0.0, 0.4), vec3(0.0, 0.3, 0.8), smoothstep(0.25, 0.5, layerPos));
        } else if (layerPos < 0.75) {
            // Upper depths - vivid teal/turquoise with higher saturation
            color = mix(vec3(0.0, 0.3, 0.8), vec3(0.0, 0.6, 1.0), smoothstep(0.5, 0.75, layerPos));
        } else {
            // Surface waters - bright cyan to create more visible layers
            color = mix(vec3(0.0, 0.6, 1.0), vec3(0.3, 0.8, 1.0), smoothstep(0.75, 1.0, layerPos));
        }

        // Add underwater caustics effect
        float causticScale = 20.0;
        float caustics = sin(depth * causticScale + time * 1.2) *
                          sin(density * causticScale + time * 1.8) * 0.5 + 0.5;
        // Stronger caustics effect
        color += vec3(0.0, 0.2, 0.3) * caustics * density * (layerPos + 0.3) * 1.5;

        // Add slight color variation based on density with stronger effect
        color *= (0.7 + density * 0.8);
    }
    // Plasma (type 7)
    else if (effectiveColorType == 7) {
        // Create noise-based distortion with higher frequency for electric look
        float noiseScale = 6.0;
        float noiseStrength = 0.3;

        // Generate base noise for large plasma patterns
        float baseNoise = fbm(vec3(depth * noiseScale,
                                 density * noiseScale * 0.4,
                                 time * 0.5));

        // Create more detailed noise for electric arcs and tendrils
        float detailNoise = noise(vec3(depth * noiseScale * 4.0,
                                     density * noiseScale * 3.0,
                                     time * 0.8));

        // Combine noise with depth for swirling effect
        float layerPos = depth + (baseNoise * 0.35 + detailNoise * 0.2) * noiseStrength;
        // Add density influence - denser regions have more electric activity
        layerPos += density * 0.2;

        // Create electric plasma color gradient
        if (layerPos < 0.25) {
            // Core plasma - deep purple to violet
            color = mix(vec3(0.3, 0.0, 0.5), vec3(0.5, 0.0, 0.8), smoothstep(0.0, 0.25, layerPos));
        } else if (layerPos < 0.5) {
            // Mid plasma - violet to magenta
            color = mix(vec3(0.5, 0.0, 0.8), vec3(0.9, 0.0, 0.9), smoothstep(0.25, 0.5, layerPos));
        } else if (layerPos < 0.75) {
            // Outer plasma - magenta to hot orange
            color = mix(vec3(0.9, 0.0, 0.9), vec3(1.0, 0.3, 0.0), smoothstep(0.5, 0.75, layerPos));
        } else {
            // Dissipating plasma - orange to yellow
            color = mix(vec3(1.0, 0.3, 0.0), vec3(1.0, 0.8, 0.0), smoothstep(0.75, 1.0, layerPos));
        }

        // Add electric arcing effect
        float arcNoise = noise(vec3(depth * 15.0, density * 15.0, time));
        float arcIntensity = pow(arcNoise, 5.0) * 2.0; // Makes arcs sharper and less frequent
        color += vec3(0.8, 0.4, 1.0) * arcIntensity * density;

        // Add pulsing glow effect
        float pulseRate = 2.0; // Speed of pulsing
        float pulseNoise = noise(vec3(depth * 2.0, density * 2.0, time * pulseRate));
        color *= 0.8 + pulseNoise * 0.4;
    }
    // Rainbow (type 8)
    else if (effectiveColorType == 8) {
        // Create noise distortion for irregular rainbow bands
        float noiseScale = 5.0;
        float noiseStrength = 0.25;

        // Base and detail noise
        float baseNoise = fbm(vec3(depth * noiseScale,
                                 density * noiseScale * 0.3,
                                 time * 0.2));
        float detailNoise = noise(vec3(depth * noiseScale * 2.0,
                                    density * noiseScale * 2.0,
                                    time * 0.2));

        // Combine with depth for flowing, irregular bands
        float bandPos = depth + (baseNoise * 0.3 + detailNoise * 0.1) * noiseStrength;
        // Add density influence
        bandPos += density * 0.1;

        // Prevent blue at the top by making sure bandPos never equals 1.0
        bandPos = min(bandPos, 0.99);

        // Create rainbow color spectrum
        if (bandPos < 0.16) {
            // Red
            color = vec3(1.0, 0.0, 0.0);
        } else if (bandPos < 0.33) {
            // Orange
            color = mix(vec3(1.0, 0.0, 0.0), vec3(1.0, 0.5, 0.0), (bandPos - 0.16) / 0.17);
        } else if (bandPos < 0.5) {
            // Yellow
            color = mix(vec3(1.0, 0.5, 0.0), vec3(1.0, 1.0, 0.0), (bandPos - 0.33) / 0.17);
        } else if (bandPos < 0.66) {
            // Green
            color = mix(vec3(1.0, 1.0, 0.0), vec3(0.0, 0.8, 0.0), (bandPos - 0.5) / 0.16);
        } else if (bandPos < 0.83) {
            // Blue
            color = mix(vec3(0.0, 0.8, 0.0), vec3(0.0, 0.0, 1.0), (bandPos - 0.66) / 0.17);
        } else {
            // Violet
            color = mix(vec3(0.0, 0.0, 1.0), vec3(0.6, 0.0, 0.8), (bandPos - 0.83) / 0.17);
        }

        // Add shimmering effect
        float shimmer = noise(vec3(depth * 20.0, density * 20.0, time * 3.0));
        color += vec3(shimmer) * 0.1;

        // Adjust brightness based on density
        color *= (0.7 + density * 0.5);
    }
    // Aurora (type 9)
    else if (effectiveColorType == 9) {
        // Create time-varying noise for aurora movement
        float timeScale = time * 0.3;
        float noiseScale = 4.0;

        // Create flowing noise patterns
        float flowNoise = fbm(vec3(depth * noiseScale,
                                 density * noiseScale,
                                 timeScale));
        float waveNoise = sin(depth * 10.0 + timeScale) *
                         cos(density * 8.0 + timeScale * 1.2) * 0.5 + 0.5;

        // Create vertical bands with distortion
        float bandPos = depth + (flowNoise * 0.3 + waveNoise * 0.1) * 0.4;
        // Density affects the band position
        bandPos += density * 0.1;

        // Auroral colors - greens, teals, purples
        if (bandPos < 0.3) {
            // Deep purple to teal
            color = mix(vec3(0.2, 0.0, 0.4), vec3(0.0, 0.4, 0.4), smoothstep(0.0, 0.3, bandPos));
        } else if (bandPos < 0.6) {
            // Teal to bright green
            color = mix(vec3(0.0, 0.4, 0.4), vec3(0.0, 1.0, 0.3), smoothstep(0.3, 0.6, bandPos));
        } else if (bandPos < 0.8) {
            // Green to light cyan
            color = mix(vec3(0.0, 1.0, 0.3), vec3(0.4, 1.0, 0.8), smoothstep(0.6, 0.8, bandPos));
        } else {
            // Cyan to purple for high regions
            color = mix(vec3(0.4, 1.0, 0.8), vec3(0.6, 0.3, 1.0), smoothstep(0.8, 1.0, bandPos));
        }

        // Add dancing light rays effect
        float rayIntensity = pow(sin(depth * 5.0 + timeScale * 2.0) *
                               sin(density * 5.0 + timeScale * 1.5), 2.0);
        color += vec3(0.2, 0.5, 0.2) * rayIntensity * (1.0 - depth) * 0.5;

        // Add subtle starfield effect for high regions
        if (depth > 0.7) {
            float stars = step(0.98, noise(vec3(depth * 50.0, density * 50.0, time * 0.5)));
            color += vec3(stars) * 0.5 * (depth - 0.7) / 0.3;
        }

        // Final color with brightness adjustment
        color *= (0.8 + density * 0.5);
    }
    // Lava (type 10)
    else if (effectiveColorType == 10) {
        // Create noise for uneven lava flow
        float noiseScale = 3.0;
        float noiseStrength = 0.3;

        // Generate churning lava patterns
        float baseNoise = fbm(vec3(depth * noiseScale,
                                 density * noiseScale * 0.2,
                                 time * 0.2));
        float detailNoise = noise(vec3(depth * noiseScale * 2.0,
                                    density * noiseScale * 1.5,
                                    time * 0.2));

        // Create bubbling, flowing lava
        float lavaPos = depth + (baseNoise * 0.4 + detailNoise * 0.2) * noiseStrength;
        // Denser regions are hotter
        lavaPos += density * 0.2;

        // Create molten lava color gradient
        if (lavaPos < 0.3) {
            // Molten core - bright orange-yellow
            color = mix(vec3(1.0, 0.8, 0.0), vec3(1.0, 0.5, 0.0), smoothstep(0.0, 0.3, lavaPos));
        } else if (lavaPos < 0.6) {
            // Mid lava - orange to red
            color = mix(vec3(1.0, 0.5, 0.0), vec3(0.9, 0.1, 0.0), smoothstep(0.3, 0.6, lavaPos));
        } else if (lavaPos < 0.85) {
            // Cooling lava - red to dark red
            color = mix(vec3(0.9, 0.1, 0.0), vec3(0.5, 0.0, 0.0), smoothstep(0.6, 0.85, lavaPos));
        } else {
            // Solidified crust - dark red to black
            color = mix(vec3(0.5, 0.0, 0.0), vec3(0.1, 0.0, 0.0), smoothstep(0.85, 1.0, lavaPos));
        }

        // Add bubbling/glowing effect using simplified approach
        float bubbleNoise = noise(vec3(depth * 15.0, density * 10.0, time * 0.5));
        float bubbleIntensity = bubbleNoise * bubbleNoise;
        color += vec3(1.0, 0.5, 0.0) * bubbleIntensity * density * 0.5;

        // Add simple cracks
        if (lavaPos > 0.7) {
            float crackNoise = noise(vec3(depth * 20.0, density * 5.0, time * 0.2));
            if (crackNoise > 0.7) {
                color += vec3(1.0, 0.4, 0.0) * (lavaPos - 0.7) * 3.0;
            }
        }

        // Add emissive glow based on temperature (higher near bottom)
        color += vec3(0.5, 0.1, 0.0) * density * pow(1.0 - lavaPos, 2.0) * 0.3;
    }

    return color;
}

// Get alpha value based on color type and render mode
float getAlphaValue(int effectiveColorType, float density, float depth, int renderMode) {
    float alphaValue;

    if (effectiveColorType == 0) {
        // Default water alpha
        if (renderMode == 1) {
            alphaValue = density * 0.5 * (1.0 + 0.6 * depth);
        } else {
            alphaValue = density * 0.18 * (1.0 + 0.5 * depth);
        }
    }
    else if (effectiveColorType >= 5 && effectiveColorType <= 7) {
        // Fire, Ocean, Plasma - higher opacity
        alphaValue = density * 0.28;
    }
    else if (effectiveColorType >= 8 && effectiveColorType <= 10) {
        // Rainbow, Aurora, Lava - highest opacity
        alphaValue = density * 0.32;
    }
    else {
        // Other types
        alphaValue = density * 0.22;
    }

    // Boost alpha in shell rendering mode
    if (renderMode == 1) {
        alphaValue *= 1.5;
    }

    return alphaValue;
}

// Apply edge highlights based on color type
vec3 applyEdgeHighlights(vec3 color, int effectiveColorType, float edge) {
    if (edge > 0.1) {
        vec3 edgeColor;
        float edgeStrength = 0.6;

        if (effectiveColorType == 5 || effectiveColorType == 10) { // Fire/Lava
            edgeColor = vec3(1.0, 0.7, 0.2); // Bright orange
            edgeStrength = 0.7;
        }
        else if (effectiveColorType == 6) { // Ocean
            edgeColor = vec3(0.7, 0.9, 1.0); // Bright cyan
        }
        else if (effectiveColorType == 7) { // Plasma
            edgeColor = vec3(0.9, 0.6, 1.0); // Bright pink
            edgeStrength = 0.7;
        }
        else if (effectiveColorType == 8) { // Rainbow
            edgeColor = vec3(1.0); // White
        }
        else if (effectiveColorType == 9) { // Aurora
            edgeColor = vec3(0.6, 1.0, 0.8); // Bright green
        }
        else {
            edgeColor = vec3(1.0); // Default white
        }

        return mix(color, edgeColor, min(edge * edgeStrength, 0.5));
    }

    return color;
}

// Add shell rendering highlights
vec3 addShellHighlights(vec3 color, int effectiveColorType) {
    if (effectiveColorType == 5) { // Fire
        return mix(color, vec3(1.0, 0.7, 0.3), 0.15);
    }
    else if (effectiveColorType == 6) { // Ocean
        return mix(color, vec3(0.4, 0.7, 1.0), 0.2);
    }
    else if (effectiveColorType == 7) { // Plasma
        return mix(color, vec3(0.8, 0.3, 1.0), 0.2);
    }
    else if (effectiveColorType == 9) { // Aurora
        return mix(color, vec3(0.0, 1.0, 0.5), 0.15);
    }
    else if (effectiveColorType == 10) { // Lava
        return mix(color, vec3(1.0, 0.3, 0.0), 0.2);
    }

    return color;
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
    float stepSize = 0.4 / float(size);  // Reduced from 0.5 for more detailed sampling
    int maxSteps = int(rayLength / stepSize) + 1;

    // Read the enter time and enter position
    float t = tEnter + 0.5 * stepSize;
    vec3 pos = rayOrigin + rayDir * t;

    // For edge highlighting
    float prevDensity = 0.0;
    vec3 prevPos = pos;

    // For water refractions
    vec3 lightDir = normalize(vec3(0.5, 1.0, 0.5));

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

        // Get cell-specific color type (if any)
        float cellColorType = texture(colorTex, texCoord).r;
        int localColorType = int(cellColorType + 0.5); // Round to nearest integer

        // Set color and alpha based on fluid style and color mapping
        vec3 color;
        float alphaValue;
        float depth = texCoord.y;  // Use y-coordinate for depth-based effects


        int effectiveColorType;
        if (localColorType > 0) {
            // If cell has its own color, use it
            effectiveColorType = localColorType;
        } else {
            // For cells with no specific color, use global color
            effectiveColorType = colorMapType;
        }



        if (localColorType == 999) {  // Special value for obstacles
            // Use a fixed color for obstacles
            color = vec3(0.4, 0.4, 0.4);  // Simple gray color
            alphaValue = density * 0.8;   // More opaque
        }
        // Pineapple parts
        else if (localColorType == 1000) {  // Pineapple body (yellow)
            // Yellow color with orange tint for pineapple
            float pineappleVar = noise(vec3(vUV.x * 15.0, vUV.y * 15.0, time * 0.05)) * 0.1;
            color = vec3(0.95, 0.8 + pineappleVar, 0.2);
            alphaValue = density * 0.98;
        }
        else if (localColorType == 1001) {  // Pineapple pattern (darker orange)
            // Orange-brown for pineapple pattern/texture
            float patternVar = noise(vec3(vUV.x * 12.0, vUV.y * 12.0, time * 0.05)) * 0.05;
            color = vec3(0.8, 0.5 + patternVar, 0.1);
            alphaValue = density * 0.98;
        }
        else if (localColorType == 1002) {  // Pineapple crown/leaves (green)
            // Bright green for crown/leaves
            float leafVar = noise(vec3(vUV.x * 20.0, vUV.y * 20.0, time * 0.1)) * 0.1;
            color = vec3(0.2, 0.8 + leafVar, 0.2);
            alphaValue = density * 0.98;
        }

        // Orange parts - pixel art style
        else if (localColorType == 1010) {  // Orange fruit
            // Bright orange color from reference
            color = vec3(1.0, 0.6, 0.2);
            alphaValue = density * 0.98;
        }
        else if (localColorType == 1011) {  // Orange stem
            // Brown stem
            color = vec3(0.5, 0.3, 0.1);
            alphaValue = density * 0.98;
        }
        else if (localColorType == 1012) {  // Orange leaf
            // Green leaf
            color = vec3(0.2, 0.6, 0.3);
            alphaValue = density * 0.98;
        }

        // Grape parts - enhanced colors with bloom effect
        else if (localColorType == 1020) {  // Grape body (dark purple)
            // Rich dark purple with bloom effect
            float bloom = pow(density, 2.0) * 0.3;
            color = vec3(0.35 + bloom, 0.1 + bloom * 0.1, 0.5 + bloom);

            // Add subtle highlight based on viewing angle
            float highlight = pow(dot(rayDir, vec3(0.0, 1.0, 0.0)), 8.0) * 0.3;
            color += vec3(highlight);

            alphaValue = density * 0.98;
        }
        else if (localColorType == 1021) {  // Grape stem (brown)
            // Natural woody brown
            color = vec3(0.4, 0.25, 0.1);
            // Add subtle wood texture
            float woodTex = noise(vec3(vUV.x * 12.0, vUV.y * 12.0, 0.0)) * 0.1;
            color *= (0.95 + woodTex);
            alphaValue = density * 0.98;
        }
        else if (localColorType == 1022) {  // Grape leaf (green)
            // Vibrant green with subtle variation
            float leafVar = noise(vec3(vUV.x * 20.0, vUV.y * 20.0, time * 0.1)) * 0.15;
            color = vec3(0.15, 0.6 + leafVar, 0.15);
            // Add vein pattern
            float vein = abs(sin(vUV.x * 20.0) * sin(vUV.y * 20.0)) * 0.2;
            color *= (1.0 - vein);
            alphaValue = density * 0.98;
        }
        else if (localColorType == 1023) {  // Grape highlight (light purple)
            // Lighter purple for highlights
            color = vec3(0.6, 0.3, 0.7);
            // Add subtle glow
            float glow = pow(density, 2.0) * 0.2;
            color += vec3(glow);
            alphaValue = density * 0.98;
        }

        // Cherry parts
        else if (localColorType == 1050) {  // Cherry fruit
            // Bright cherry red
            color = vec3(0.9, 0.1, 0.1);
            // Add subtle highlight for roundness effect
            float highlight = pow(density, 3.0) * 0.3;
            color += vec3(highlight);
            alphaValue = density * 0.98;
        }
        else if (localColorType == 1051) {  // Cherry stem
            // Brown-green stem
            color = vec3(0.25, 0.35, 0.05);
            alphaValue = density * 0.98;
        }
        else if (localColorType == 1052) {  // Cherry leaf
            // Bright green leaf
            float leafVar = noise(vec3(vUV.x * 20.0, vUV.y * 20.0, time * 0.1)) * 0.15;
            color = vec3(0.3, 0.8 + leafVar, 0.1);
            alphaValue = density * 0.98;
        }
        // Pear parts
        else if (localColorType == 1040) {  // Pear fruit
            // Yellow-green pear color with subtle variation
            float pearVar = noise(vec3(vUV.x * 15.0, vUV.y * 15.0, time * 0.05)) * 0.1;
            color = vec3(0.8, 0.9 + pearVar, 0.3 + pearVar);

            // Add subtle highlight
            float highlight = pow(density, 3.0) * 0.2;
            color += vec3(highlight);

            // Add subtle texture pattern
            float texture = noise(vec3(vUV.x * 25.0, vUV.y * 25.0, time * 0.01)) * 0.05;
            color.g += texture;

            alphaValue = density * 0.98;
        }
        else if (localColorType == 1041) {  // Pear stem
            // Brown stem with texture
            color = vec3(0.4, 0.25, 0.1);
            float stemTex = noise(vec3(vUV.x * 20.0, vUV.y * 5.0, time * 0.05)) * 0.1;
            color *= (0.95 + stemTex);

            alphaValue = density * 0.98;
        }
        else {
            // Get base color for this fluid type
            color = getFluidColor(effectiveColorType, depth, density, edge);

            // Get alpha value
            alphaValue = getAlphaValue(effectiveColorType, density, depth, renderMode);

            // Apply edge highlights
            color = applyEdgeHighlights(color, effectiveColorType, edge);

            // Add special shell rendering effects
            if (renderMode == 1) {
                color = addShellHighlights(color, effectiveColorType);
            }

            // Add Fresnel effect only for default water
            if (effectiveColorType == 0) {
                float fresnelEffect = fresnel(abs(dot(rayDir, vec3(0.0, 1.0, 0.0))), 0.05);
                color = mix(color, vec3(1.0), fresnelEffect * 0.4);
            }

            // Special top surface handling to prevent blue artifacts
            if (texCoord.y > 0.95 && effectiveColorType > 0) {
                // Near the top surface, gradually fade to the scheme's own color instead of default blue
                float topFactor = (texCoord.y - 0.95) / 0.05;
                // Reduce any blue component to prevent blue artifacts
                color.b = mix(color.b, min(color.b, 0.5), topFactor);
            }
        }


        // Alpha Blending
        finalColor.rgb += (1.0 - finalColor.a) * alphaValue * color;
        finalColor.a += (1.0 - finalColor.a) * alphaValue;

        // Exit earlier if the alpha is already large enough
        if (finalColor.a > 0.80)
            break;

        t += stepSize;
        pos += rayDir * stepSize;
    }

    // Post-processing effects
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
