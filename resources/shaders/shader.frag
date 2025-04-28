#version 330 core
// out vec4 fragColor;

// // These are for the original FEM project
// in vec4 normal_worldSpace;
// in vec4 position_worldSpace;

// uniform int wire = 0;
// uniform float red = 1.0;
// uniform float green = 1.0;
// uniform float blue = 1.0;
// uniform float alpha = 1.0;

// // These are for our final project
// in vec2 vUV;

// uniform float gray; // For Voxel
// uniform sampler3D densityTex; // For RayMarching
// uniform sampler3D colorTex;   // For per-cell color mapping
// uniform mat4 invViewProj;
// uniform float time; // For time-based animation effects

// // Fluid style: 0=gray, 1=blue water
// uniform int fluidStyle = 1;  // Default to blue water

// // Color map types - extended with new gradient effects
// uniform int colorMapType = 0; // 0: default, 1: blue, 2: purple, 3: cyan-yellow, 4: orange-grey
//                              // 5: fire gradient, 6: ocean gradient, 7: plasma gradient
//                              // 8: rainbow, 9: aurora, 10: lava

// uniform int renderMode = 0; // 0: volume, 1: shell

// // Define the position of the box - increased from 0.5 to 2.0 for larger display
// const vec3 boxMin = vec3(-2.0);
// const vec3 boxMax = vec3(2.0);

// // Set the step parameters in ray marching
// // If the simulation is not detailed enough, we should decrease stepSize or increase maxSteps
// uniform int size;
// float voxelLength = 1.0 / float(size);

// // Simple hash function for noise generation
// float hash(float n) {
//     return fract(sin(n) * 43758.5453);
// }

// // 3D noise function for natural randomness
// float noise(vec3 x) {
//     vec3 p = floor(x);
//     vec3 f = fract(x);
//     f = f * f * (3.0 - 2.0 * f);

//     float n = p.x + p.y * 157.0 + 113.0 * p.z;
//     return mix(mix(mix(hash(n + 0.0), hash(n + 1.0), f.x),
//                    mix(hash(n + 157.0), hash(n + 158.0), f.x), f.y),
//                mix(mix(hash(n + 113.0), hash(n + 114.0), f.x),
//                    mix(hash(n + 270.0), hash(n + 271.0), f.x), f.y), f.z);
// }

// // Improved Fractal Brownian Motion for more natural noise with stronger turbulence
// float fbm(vec3 x, float turbulence) {
//     float v = 0.0;
//     float a = 0.5;
//     vec3 shift = vec3(100);

//     // Rotate to reduce axial bias
//     mat3 rot = mat3(cos(0.5), sin(0.5), 0,
//                     -sin(0.5), cos(0.5), 0,
//                     0, 0, 1);

//     for (int i = 0; i < 6; ++i) { // More octaves for greater detail
//         v += a * noise(x);
//         x = rot * x * 2.0 + shift;
//         a *= 0.5;
//     }
//     // Apply turbulence factor to control noise strength
//     return v * turbulence;
// }

// // Turbulent noise with more dramatic variation
// float turbulence(vec3 p, float freq, float lacunarity, float gain) {
//     float sum = 0.0;
//     float amp = 1.0;
//     vec3 pp = p * freq;

//     for(int i = 0; i < 6; i++) {
//         sum += abs(noise(pp) * amp);
//         amp *= gain;
//         pp *= lacunarity;
//     }
//     return sum;
// }

// // Curl noise for more natural flow effects
// vec3 curlNoise(vec3 p, float time) {
//     // Small delta for numerical derivatives
//     float epsilon = 0.001;

//     // Calculate vorticity through curl of a noise field
//     float n1, n2, n3, n4, n5, n6;

//     // Add time influence for animation
//     vec3 pTime = p + vec3(0.0, time * 0.1, 0.0);

//     // Derivative X
//     n1 = noise(pTime + vec3(0, epsilon, 0));
//     n2 = noise(pTime - vec3(0, epsilon, 0));
//     float dydx = (n1 - n2) / (2.0 * epsilon);

//     n1 = noise(pTime + vec3(0, 0, epsilon));
//     n2 = noise(pTime - vec3(0, 0, epsilon));
//     float dzdx = (n1 - n2) / (2.0 * epsilon);

//     // Derivative Y
//     n3 = noise(pTime + vec3(epsilon, 0, 0));
//     n4 = noise(pTime - vec3(epsilon, 0, 0));
//     float dxdy = (n3 - n4) / (2.0 * epsilon);

//     n3 = noise(pTime + vec3(0, 0, epsilon));
//     n4 = noise(pTime - vec3(0, 0, epsilon));
//     float dzdy = (n3 - n4) / (2.0 * epsilon);

//     // Derivative Z
//     n5 = noise(pTime + vec3(epsilon, 0, 0));
//     n6 = noise(pTime - vec3(epsilon, 0, 0));
//     float dxdz = (n5 - n6) / (2.0 * epsilon);

//     n5 = noise(pTime + vec3(0, epsilon, 0));
//     n6 = noise(pTime - vec3(0, epsilon, 0));
//     float dydz = (n5 - n6) / (2.0 * epsilon);

//     return vec3(dydz - dzdy, dzdx - dxdz, dxdy - dydx);
// }

// // Using AABB to detect the intersection
// bool intersectBox(vec3 origin, vec3 dir, out float tMin, out float tMax) {
//     vec3 invDir = 1.0 / dir;
//     vec3 t0s = (boxMin - origin) * invDir;
//     vec3 t1s = (boxMax - origin) * invDir;

//     vec3 tSmalls = min(t0s, t1s);
//     vec3 tBigs   = max(t0s, t1s);

//     tMin = max(max(tSmalls.x, tSmalls.y), tSmalls.z);
//     tMax = min(min(tBigs.x, tBigs.y), tBigs.z);

//     return tMax >= tMin;
// }

// // Helper functions for water effects
// float fresnel(float cosTheta, float F0) {
//     return F0 + (1.0 - F0) * pow(1.0 - cosTheta, 5.0);
// }

// // Create swirly flow distortion field
// vec3 distortPosition(vec3 position, int colorType, float time) {
//     // Base parameters
//     float distortScale = 3.0;
//     float timeScale = time * 0.3;

//     // Customize distortion by color type
//     if (colorType == 5) { // Fire
//         // Rising hot currents for fire
//         vec3 curl = curlNoise(position * distortScale, timeScale);
//         return position + vec3(curl.x * 0.05, curl.y * 0.15, curl.z * 0.05);
//     }
//     else if (colorType == 6) { // Ocean
//         // Horizontal currents for water
//         vec3 curl = curlNoise(position * distortScale * 0.7, timeScale * 0.5);
//         return position + vec3(curl.x * 0.08, curl.y * 0.03, curl.z * 0.08);
//     }
//     else if (colorType == 7) { // Plasma
//         // Chaotic flow for plasma
//         vec3 curl = curlNoise(position * distortScale * 1.5, timeScale * 1.5);
//         return position + curl * 0.08;
//     }
//     else if (colorType == 8) { // Rainbow
//         // Gentle turbulence for rainbow
//         vec3 curl = curlNoise(position * distortScale * 0.6, timeScale * 0.4);
//         return position + curl * 0.04;
//     }
//     else if (colorType == 9) { // Aurora
//         // Flowing bands for aurora
//         vec3 curl = curlNoise(position * vec3(1.5, 0.5, 1.5), timeScale * 0.7);
//         return position + vec3(curl.x * 0.1, curl.y * 0.02, curl.z * 0.1);
//     }
//     else if (colorType == 10) { // Lava
//         // Slow bubbling for lava
//         vec3 curl = curlNoise(position * distortScale * 0.8, timeScale * 0.3);
//         return position + curl * 0.07;
//     }
//     else {
//         // Default mild distortion
//         vec3 curl = curlNoise(position * distortScale, timeScale);
//         return position + curl * 0.03;
//     }
// }

// // Create non-linear layer position based on coordinates and fluid type
// float createLayerPosition(vec3 texCoord, float density, int colorType, float time) {
//     // Base parameters to be customized per color type
//     float noiseScale = 4.0;
//     float noiseStrength = 0.4; // Increased for more dramatic distortion
//     float verticalInfluence = 0.5;
//     float densityInfluence = 0.2;

//     // Adjust parameters by color type
//     if (colorType == 5) { // Fire
//         noiseScale = 5.0;
//         noiseStrength = 0.45;
//         verticalInfluence = 0.4;
//         densityInfluence = 0.3;
//     }
//     else if (colorType == 6) { // Ocean
//         noiseScale = 3.0;
//         noiseStrength = 0.5;
//         verticalInfluence = 0.6;
//         densityInfluence = -0.15; // Negative: higher density appears deeper
//     }
//     else if (colorType == 7) { // Plasma
//         noiseScale = 6.0;
//         noiseStrength = 0.4;
//         verticalInfluence = 0.5;
//         densityInfluence = 0.25;
//     }
//     else if (colorType == 8) { // Rainbow
//         noiseScale = 4.0;
//         noiseStrength = 0.35;
//         verticalInfluence = 0.6;
//         densityInfluence = 0.1;
//     }
//     else if (colorType == 9) { // Aurora
//         noiseScale = 3.5;
//         noiseStrength = 0.45;
//         verticalInfluence = 0.6;
//         densityInfluence = 0.15;
//     }
//     else if (colorType == 10) { // Lava
//         noiseScale = 3.0;
//         noiseStrength = 0.5;
//         verticalInfluence = 0.4;
//         densityInfluence = 0.25;
//     }

//     // Create time-varying turbulence
//     float timeMod = sin(time * 0.2) * 0.5 + 0.5; // Oscillating turbulence
//     float turbFactor = 1.0 + timeMod * 0.2;

//     // Generate noise using position and time
//     vec3 noisePos = texCoord * noiseScale;
//     noisePos.y *= 0.3; // Less variation vertically

//     // Create three different noise components at different scales
//     float baseNoise = fbm(noisePos, turbFactor);
//     float detailNoise = noise(vec3(texCoord.x * noiseScale * 2.0,
//                                   texCoord.y * noiseScale * 2.0,
//                                   texCoord.z * noiseScale * 2.0 + time * 0.3));

//     // Add extra turbulence for more natural, chaotic look
//     float turbNoise = turbulence(noisePos + vec3(0, time * 0.1, 0), 1.0, 2.0, 0.5) * 0.1;

//     // Combine with vertical position for basic layering
//     float layerPos = texCoord.y * verticalInfluence +
//                     (baseNoise * 0.4 + detailNoise * 0.2 + turbNoise) * noiseStrength;

//     // Add density influence - varies by type
//     layerPos += density * densityInfluence;

//     // Add subtle vortex/swirl at high density areas
//     if (density > 0.5) {
//         float vortexNoise = noise(vec3(texCoord.xz * 8.0, time * 0.5));
//         layerPos += vortexNoise * density * 0.1;
//     }

//     // Ensure layerPos stays in 0-1 range
//     return clamp(layerPos, 0.0, 0.99);
// }

// // Color mapping function that handles all types and prevents blue lines at top
// vec3 getFluidColor(int colorType, float layerPos, float depth, float density, vec3 texCoord, float time) {
//     // Common variables
//     vec3 color;
//     bool isTopArea = depth > 0.92; // Top 8% of volume

//     // Default water
//     if (colorType == 0) {
//         vec3 waterBaseColor = vec3(0.1, 0.4, 0.8);
//         vec3 waterDeepColor = vec3(0.0, 0.1, 0.3);
//         vec3 waterHighlightColor = vec3(0.8, 0.9, 1.0);

//         float depthEffect = pow(depth, 1.2);
//         color = mix(waterDeepColor, waterBaseColor, depthEffect);

//         // Add foam effect to high density areas
//         color = mix(color, waterHighlightColor, pow(density, 2.5) * 0.6);

//         // Add surface reflection at top of fluid
//         float surfaceFactor = smoothstep(0.75, 1.0, depth);
//         color = mix(color, vec3(0.5, 0.7, 1.0), surfaceFactor * 0.6);

//         return color;
//     }

//     // Cyan-blue (type 1)
//     else if (colorType == 1) {
//         color = mix(vec3(0.0, 0.9, 1.0), vec3(0.0, 0.1, 0.6), pow(1.0-layerPos, 1.5));

//         // Top area handling
//         if (isTopArea) {
//             color = mix(color, vec3(0.0, 0.7, 0.9), (depth - 0.92) / 0.08);
//         }
//     }

//     // Purple electric (type 2)
//     else if (colorType == 2) {
//         color = mix(vec3(0.7, 0.0, 1.0), vec3(0.1, 0.0, 0.4), pow(1.0-layerPos, 1.3));

//         // Top area handling
//         if (isTopArea) {
//             color = mix(color, vec3(0.5, 0.0, 0.7), (depth - 0.92) / 0.08);
//             color.b *= 0.5; // Reduce blue
//         }
//     }

//     // Cyan-yellow (type 3)
//     else if (colorType == 3) {
//         color = mix(vec3(0.0, 0.9, 0.9), vec3(1.0, 0.9, 0.1), pow(layerPos, 1.2));

//         // Top area handling
//         if (isTopArea) {
//             color = mix(color, vec3(0.9, 0.8, 0.1), (depth - 0.92) / 0.08);
//             color.b *= 0.3; // Reduce blue significantly
//         }
//     }

//     // Orange-grey (type 4)
//     else if (colorType == 4) {
//         color = mix(vec3(1.0, 0.4, 0.0), vec3(0.6, 0.6, 0.6), pow(layerPos, 1.1));

//         // Top area handling
//         if (isTopArea) {
//             color = mix(color, vec3(0.7, 0.5, 0.3), (depth - 0.92) / 0.08);
//             color.b *= 0.5; // Reduce blue
//         }
//     }

//     // Fire (type 5)
//     else if (colorType == 5) {
//         // Create multiple color bands with noise-influenced transitions
//         if (layerPos < 0.25) {
//             // Inner core - bright yellow-white
//             color = mix(vec3(1.0, 0.9, 0.5), vec3(1.0, 0.7, 0.0), smoothstep(0.0, 0.25, layerPos));
//         } else if (layerPos < 0.5) {
//             // Mid flame - yellow to orange
//             color = mix(vec3(1.0, 0.7, 0.0), vec3(1.0, 0.4, 0.0), smoothstep(0.25, 0.5, layerPos));
//         } else if (layerPos < 0.75) {
//             // Outer flame - orange to red
//             color = mix(vec3(1.0, 0.4, 0.0), vec3(0.9, 0.1, 0.0), smoothstep(0.5, 0.75, layerPos));
//         } else {
//             // Fading edges - red to dark red
//             color = mix(vec3(0.9, 0.1, 0.0), vec3(0.4, 0.0, 0.0), smoothstep(0.75, 1.0, layerPos));
//         }

//         // Add flickering effect
//         float flicker = noise(vec3(texCoord.x * 10.0, texCoord.z * 10.0, time * 2.0));
//         color *= 0.8 + flicker * 0.4;

//         // Add turbulent glow near dense areas
//         float glowNoise = noise(texCoord * 15.0 + vec3(0, time * 0.5, 0));
//         color += vec3(0.3, 0.1, 0.0) * glowNoise * density * 0.6;

//         // Top area handling - ensure no blue
//         if (isTopArea) {
//             color = mix(color, vec3(0.7, 0.1, 0.0), (depth - 0.92) / 0.08);
//             color.b = 0.0; // No blue at all
//         }
//     }

//     // Ocean (type 6)
//     else if (colorType == 6) {
//         // Create more dramatic ocean bands with enhanced contrast
//         if (layerPos < 0.2) {
//             // Deep ocean - dark blue
//             color = vec3(0.0, 0.0, 0.2);
//         } else if (layerPos < 0.4) {
//             // Mid-depths - rich blue
//             color = mix(vec3(0.0, 0.0, 0.3), vec3(0.0, 0.2, 0.6), smoothstep(0.2, 0.4, layerPos));
//         } else if (layerPos < 0.65) {
//             // Upper-mid depths - blue to teal
//             color = mix(vec3(0.0, 0.2, 0.6), vec3(0.0, 0.5, 0.8), smoothstep(0.4, 0.65, layerPos));
//         } else if (layerPos < 0.85) {
//             // Upper depths - teal to cyan
//             color = mix(vec3(0.0, 0.5, 0.8), vec3(0.0, 0.7, 0.9), smoothstep(0.65, 0.85, layerPos));
//         } else {
//             // Surface - cyan to light blue-green
//             color = mix(vec3(0.0, 0.7, 0.9), vec3(0.2, 0.8, 0.8), smoothstep(0.85, 1.0, layerPos));
//         }

//         // Add underwater caustics effect
//         float causticScale = 20.0;
//         float caustics = sin(texCoord.x * causticScale + time * 0.7) *
//                         sin(texCoord.z * causticScale + time * 1.0) * 0.5 + 0.5;
//         // Apply caustics more strongly in mid-layers
//         float causticIntensity = (1.0 - abs(layerPos - 0.5) * 2.0) * 0.6;
//         color += vec3(0.0, 0.2, 0.3) * caustics * density * causticIntensity;

//         // Add waving surface effect at top
//         if (depth > 0.8) {
//             float waveEffect = sin(texCoord.x * 10.0 + time) * cos(texCoord.z * 8.0 + time * 1.2) * 0.5 + 0.5;
//             color += vec3(0.1, 0.2, 0.2) * waveEffect * (depth - 0.8) / 0.2;
//         }

//         // Top area handling - ensure cyan rather than blue
//         if (isTopArea) {
//             color = mix(color, vec3(0.0, 0.7, 0.7), (depth - 0.92) / 0.08);
//         }
//     }

//     // Plasma (type 7)
//     else if (colorType == 7) {
//         // Create electric plasma gradient with more bands
//         if (layerPos < 0.2) {
//             // Core plasma - deep purple
//             color = vec3(0.3, 0.0, 0.5);
//         } else if (layerPos < 0.4) {
//             // Inner plasma - purple to violet
//             color = mix(vec3(0.3, 0.0, 0.5), vec3(0.5, 0.0, 0.8), smoothstep(0.2, 0.4, layerPos));
//         } else if (layerPos < 0.6) {
//             // Mid plasma - violet to magenta
//             color = mix(vec3(0.5, 0.0, 0.8), vec3(0.9, 0.0, 0.8), smoothstep(0.4, 0.6, layerPos));
//         } else if (layerPos < 0.8) {
//             // Outer plasma - magenta to hot pink/orange
//             color = mix(vec3(0.9, 0.0, 0.8), vec3(1.0, 0.2, 0.2), smoothstep(0.6, 0.8, layerPos));
//         } else {
//             // Edge plasma - orange to yellow
//             color = mix(vec3(1.0, 0.2, 0.2), vec3(1.0, 0.7, 0.0), smoothstep(0.8, 1.0, layerPos));
//         }

//         // Add electric arcing effect
//         float arcNoise = noise(vec3(texCoord.x * 15.0, texCoord.y * 15.0, texCoord.z * 15.0 + time));
//         float arcIntensity = pow(arcNoise, 5.0) * 2.0; // Make arcs sharper and less frequent
//         color += vec3(0.8, 0.4, 1.0) * arcIntensity * density;

//         // Add pulsing glow effect
//         float pulseRate = 2.0; // Speed of pulsing
//         float pulseNoise = noise(vec3(texCoord.x * 2.0, texCoord.y * 2.0, texCoord.z * 2.0 + time * pulseRate));
//         color *= 0.8 + pulseNoise * 0.4;

//         // Top area handling - reduce blue
//         if (isTopArea) {
//             color = mix(color, vec3(0.9, 0.2, 0.5), (depth - 0.92) / 0.08);
//             color.b *= 0.4; // Reduce blue
//         }
//     }

//     // Rainbow (type 8)
//     else if (colorType == 8) {
//         // Create smoother rainbow spectrum with 7 bands
//         float band = layerPos * 7.0; // 7 color bands
//         int bandIndex = int(band);
//         float bandMix = fract(band);

//         // Define rainbow colors
//         vec3 rainbowColors[8];
//         rainbowColors[0] = vec3(0.8, 0.0, 0.0); // Red
//         rainbowColors[1] = vec3(1.0, 0.4, 0.0); // Orange
//         rainbowColors[2] = vec3(1.0, 0.8, 0.0); // Yellow
//         rainbowColors[3] = vec3(0.0, 0.8, 0.0); // Green
//         rainbowColors[4] = vec3(0.0, 0.4, 0.8); // Blue
//         rainbowColors[5] = vec3(0.2, 0.0, 0.8); // Indigo
//         rainbowColors[6] = vec3(0.6, 0.0, 0.8); // Violet
//         rainbowColors[7] = vec3(0.8, 0.0, 0.4); // Return to reddish-purple

//         // Mix between appropriate bands
//         color = mix(rainbowColors[bandIndex], rainbowColors[bandIndex+1], bandMix);

//         // Add shimmer effect
//         float shimmer = noise(vec3(texCoord.x * 20.0, texCoord.z * 20.0, time * 2.0));
//         color += vec3(shimmer) * 0.1;

//         // Top area handling - shift towards red/violet
//         if (isTopArea) {
//             color = mix(color, vec3(0.7, 0.0, 0.4), (depth - 0.92) / 0.08);
//         }
//     }

//     // Aurora (type 9)
//     else if (colorType == 9) {
//         // Create time-varying noise for aurora movement
//         float timeScale = time * 0.3;

//         // Layer position influences base color
//         if (layerPos < 0.25) {
//             // Deep base - purple to teal
//             color = mix(vec3(0.1, 0.0, 0.3), vec3(0.0, 0.3, 0.3), smoothstep(0.0, 0.25, layerPos));
//         } else if (layerPos < 0.5) {
//             // Mid aurora - teal to green
//             color = mix(vec3(0.0, 0.3, 0.3), vec3(0.0, 0.7, 0.2), smoothstep(0.25, 0.5, layerPos));
//         } else if (layerPos < 0.75) {
//             // Upper aurora - green to light cyan
//             color = mix(vec3(0.0, 0.7, 0.2), vec3(0.3, 0.8, 0.7), smoothstep(0.5, 0.75, layerPos));
//         } else {
//             // Top aurora - light cyan to purple
//             color = mix(vec3(0.3, 0.8, 0.7), vec3(0.5, 0.2, 0.8), smoothstep(0.75, 1.0, layerPos));
//         }

//         // Add dancing light rays effect
//         float rayX = sin(texCoord.x * 5.0 + timeScale * 2.0);
//         float rayZ = sin(texCoord.z * 5.0 + timeScale * 1.5);
//         float rayIntensity = pow(rayX * rayZ, 2.0) * (1.0 - texCoord.y);

//         // Get ray color based on position - different rays have different colors
//         vec3 rayColor;
//         float rayPos = noise(vec3(texCoord.x * 3.0, 0.0, texCoord.z * 3.0));
//         if (rayPos < 0.33) {
//             rayColor = vec3(0.0, 0.8, 0.2); // Green ray
//         } else if (rayPos < 0.66) {
//             rayColor = vec3(0.0, 0.5, 0.8); // Blue ray
//         } else {
//             rayColor = vec3(0.5, 0.0, 0.8); // Purple ray
//         }

//         color += rayColor * rayIntensity * 0.5;

//         // Add subtle starfield effect for high regions
//         if (depth > 0.7) {
//             float stars = step(0.98, noise(vec3(texCoord.x * 50.0, texCoord.y * 50.0, texCoord.z * 50.0)));
//             color += vec3(stars) * 0.5 * (depth - 0.7) / 0.3;
//         }

//         // Top area handling - shift towards green with minimal blue
//         if (isTopArea) {
//             color = mix(color, vec3(0.2, 0.8, 0.2), (depth - 0.92) / 0.08);
//             color.b *= 0.5; // Reduce blue
//         }
//     }

//     // Lava (type 10)
//     else if (colorType == 10) {
//         // Create bubbling, flowing lava with more detailed bands
//         if (layerPos < 0.2) {
//             // Molten core - bright orange-yellow
//             color = mix(vec3(1.0, 0.9, 0.0), vec3(1.0, 0.7, 0.0), smoothstep(0.0, 0.2, layerPos));
//         } else if (layerPos < 0.4) {
//             // Inner lava - orange
//             color = mix(vec3(1.0, 0.7, 0.0), vec3(1.0, 0.5, 0.0), smoothstep(0.2, 0.4, layerPos));
//         } else if (layerPos < 0.6) {
//             // Mid lava - orange to red
//             color = mix(vec3(1.0, 0.5, 0.0), vec3(0.9, 0.2, 0.0), smoothstep(0.4, 0.6, layerPos));
//         } else if (layerPos < 0.8) {
//             // Cooling lava - red to dark red
//             color = mix(vec3(0.9, 0.2, 0.0), vec3(0.6, 0.0, 0.0), smoothstep(0.6, 0.8, layerPos));
//         } else {
//             // Crust - dark red to black
//             color = mix(vec3(0.6, 0.0, 0.0), vec3(0.2, 0.0, 0.0), smoothstep(0.8, 1.0, layerPos));
//         }

//         // Add bubbling/glowing effect
//         float bubbleNoise = noise(vec3(texCoord.x * 15.0, texCoord.y * 10.0, texCoord.z * 15.0 + time * 0.3));
//         float bubbleIntensity = bubbleNoise * bubbleNoise; // Squaring gives more contrast

//         // More bubbles in lower sections
//         float bubbleFactor = (1.0 - layerPos) * 0.7;
//         color += vec3(1.0, 0.5, 0.0) * bubbleIntensity * density * bubbleFactor;

//         // Add cracks in cooling crust
//         if (layerPos > 0.7) {
//             float crackNoise = noise(vec3(texCoord.x * 20.0, texCoord.y * 5.0, texCoord.z * 20.0));
//             if (crackNoise > 0.7) {
//                 color += vec3(1.0, 0.4, 0.0) * (layerPos - 0.7) * 3.0;
//             }
//         }

//         // Add emissive glow that decreases with depth
//         color += vec3(0.5, 0.1, 0.0) * density * pow(1.0 - layerPos, 2.0) * 0.3;

//         // Top area handling - ensure no blue
//         if (isTopArea) {
//             color = mix(color, vec3(0.3, 0.0, 0.0), (depth - 0.92) / 0.08);
//             color.b = 0.0; // No blue at all
//         }
//     }

//     return color;
// }

// void main() {
//     // Step1: vUV [0, 1] to [-1, 1] (NDC coord)
//     vec2 screenPos = vUV * 2.0 - 1.0;

//     // Step2: Construct near and far planes' points
//     vec4 pointNear = invViewProj * vec4(screenPos, 0.0, 1.0);
//     vec4 pointFar = invViewProj * vec4(screenPos, 1.0, 1.0);

//     // Step3: Construct the ray
//     vec3 rayOrigin = pointNear.xyz / pointNear.w;
//     vec3 rayTarget = pointFar.xyz / pointFar.w;
//     vec3 rayDir = normalize(rayTarget - rayOrigin);

//     // The time that the ray enters the cube and exits the cube
//     float tEnter, tExit;
//     // If not intersect, discard this ray
//     if (!intersectBox(rayOrigin, rayDir, tEnter, tExit)) discard;

//     // Ray Marching Logic starts here
//     vec4 finalColor = vec4(0.0);

//     // What we want for steps: (tMax - tMin) / stepSize ≈ maxSteps
//     float rayLength = tExit - tEnter;
//     float stepSize = 0.4 / float(size);  // Reduced from 0.5 for more detailed sampling
//     int maxSteps = int(rayLength / stepSize) + 1;

//     // Read the enter time and enter position
//     float t = tEnter + 0.5 * stepSize;
//     vec3 pos = rayOrigin + rayDir * t;

//     // For edge highlighting
//     float prevDensity = 0.0;
//     vec3 prevPos = pos;

//     // For water refractions
//     vec3 lightDir = normalize(vec3(0.5, 1.0, 0.5));
//     float NdotL = max(dot(rayDir, lightDir), 0.0);

//     for (int i = 0; i < maxSteps && t < tExit; ++i) {
//         // Map the pos to [0, 1] to get the related density
//         vec3 texCoord = (pos - boxMin) / (boxMax - boxMin);

//         // If texCoord out of boundaries, exit
//         if (any(lessThan(texCoord, vec3(0.0))) || any(greaterThan(texCoord, vec3(1.0))))
//             break;

//         // Get the density
//         float density = texture(densityTex, texCoord).r;

//         // Skip low density regions for efficiency
//         if (density < 0.01) {
//             t += stepSize * 2.0;
//             pos += rayDir * stepSize * 2.0;
//             continue;
//         }

//         // Calculate edge factor (density gradient)
//         float edge = 0.0;
//         if (i > 0) {
//             edge = abs(density - prevDensity) * 15.0;
//         }
//         prevDensity = density;
//         prevPos = pos;

//         // Get cell-specific color type (if any)
//         float cellColorType = texture(colorTex, texCoord).r;
//         int localColorType = int(cellColorType + 0.5); // Round to nearest integer

//         // Determine effective color type
//         int effectiveColorType;
//         if (localColorType > 0) {
//             // If cell has its own color, use it
//             effectiveColorType = localColorType;
//         } else if (density < 0.05) {
//             // For very low density areas, use default coloring
//             effectiveColorType = 0;
//         } else {
//             // For cells with no color but sufficient density, use global color
//             effectiveColorType = colorMapType;
//         }

//         // Set color and alpha based on fluid style and color mapping
//         vec3 color;
//         float alphaValue;
//         float depth = texCoord.y;  // Use y-coordinate for depth-based effects

//         if (fluidStyle == 0) {  // Gray (original)
//             color = vec3(1.0 - density);
//             alphaValue = density * 0.15;
//         } else if (fluidStyle == 1) {  // Color fluid with enhanced effects
//             // Distort the position based on color type for more natural flow
//             vec3 distortedPos = distortPosition(texCoord, effectiveColorType, time);

//             // Calculate layer position with multiple noise functions for more natural banding
//             float layerPos = createLayerPosition(distortedPos, density, effectiveColorType, time);

//             // Get color based on the layer position and type
//             color = getFluidColor(effectiveColorType, layerPos, depth, density, texCoord, time);

//             // Determine alpha value based on color type
//             if (effectiveColorType == 0) {
//                 // Default water alpha
//                 if (renderMode == 1) {
//                     alphaValue = density * 0.5 * (1.0 + 0.6 * depth);
//                 } else {
//                     alphaValue = density * 0.18 * (1.0 + 0.5 * depth);
//                 }
//             }
//             else if (effectiveColorType >= 5 && effectiveColorType <= 7) {
//                 // Fire, Ocean, Plasma - higher opacity
//                 alphaValue = density * 0.28;
//             }
//             else if (effectiveColorType >= 8 && effectiveColorType <= 10) {
//                 // Rainbow, Aurora, Lava - highest opacity
//                 alphaValue = density * 0.32;
//             }
//             else {
//                 // Other types
//                 alphaValue = density * 0.22;
//             }

//             // Enhanced edge highlights - customize by color type
//             if (edge > 0.1) {
//                 vec3 edgeColor;
//                 float edgeStrength = 0.6;

//                 if (effectiveColorType == 5 || effectiveColorType == 10) { // Fire/Lava
//                     edgeColor = vec3(1.0, 0.7, 0.2); // Bright orange
//                     edgeStrength = 0.7;
//                 }
//                 else if (effectiveColorType == 6) { // Ocean
//                     edgeColor = vec3(0.4, 0.9, 1.0); // Bright cyan
//                 }
//                 else if (effectiveColorType == 7) { // Plasma
//                     edgeColor = vec3(1.0, 0.5, 1.0); // Bright pink
//                     edgeStrength = 0.7;
//                 }
//                 else if (effectiveColorType == 8) { // Rainbow
//                     edgeColor = vec3(1.0); // White
//                 }
//                 else if (effectiveColorType == 9) { // Aurora
//                     edgeColor = vec3(0.4, 1.0, 0.6); // Bright green
//                 }
//                 else {
//                     edgeColor = vec3(1.0); // Default white
//                 }

//                 color = mix(color, edgeColor, min(edge * edgeStrength, 0.5));
//             }

//             // Add Fresnel effect only for default water
//             if (effectiveColorType == 0) {
//                 float fresnelEffect = fresnel(abs(dot(rayDir, vec3(0.0, 1.0, 0.0))), 0.05);
//                 color = mix(color, vec3(1.0), fresnelEffect * 0.4);
//             }

//             // Boost alpha in shell rendering mode
//             if (renderMode == 1) {
//                 alphaValue *= 1.5;

//                 // Add extra highlighting in shell mode
//                 if (effectiveColorType == 5) { // Fire
//                     color = mix(color, vec3(1.0, 0.7, 0.3), 0.15);
//                 }
//                 else if (effectiveColorType == 6) { // Ocean
//                     color = mix(color, vec3(0.3, 0.7, 1.0), 0.15);
//                 }
//                 else if (effectiveColorType == 7) { // Plasma
//                     color = mix(color, vec3(0.8, 0.3, 1.0), 0.15);
//                 }
//                 else if (effectiveColorType == 9) { // Aurora
//                     color = mix(color, vec3(0.0, 1.0, 0.5), 0.15);
//                 }
//                 else if (effectiveColorType == 10) { // Lava
//                     color = mix(color, vec3(1.0, 0.3, 0.0), 0.15);
//                 }
//             }
//         }

//         // Alpha Blending
//         finalColor.rgb += (1.0 - finalColor.a) * alphaValue * color;
//         finalColor.a += (1.0 - finalColor.a) * alphaValue;

//         // Exit earlier if the alpha is already large enough
//         if (finalColor.a > 0.95)
//             break;

//         t += stepSize;
//         pos += rayDir * stepSize;
//     }

//     // Post-processing effects based on fluid style
//     if (finalColor.a > 0.0) {
//         if (fluidStyle == 1) {  // Blue water
//             // Add subtle tint to the entire scene - but ONLY for default color type
//             if (colorMapType == 0) {
//                 finalColor.rgb = mix(finalColor.rgb, vec3(0.4, 0.6, 0.9), 0.05);
//             }

//             // Add subtle caustics effect
//             float caustic = sin(vUV.x * 20.0) * sin(vUV.y * 20.0) * 0.5 + 0.5;
//             finalColor.rgb += vec3(0.0, 0.0, 0.2) * caustic * 0.05;
//         }

//         // Add subtle vignette effect for all styles
//         float vignette = 1.0 - dot(vUV - 0.5, vUV - 0.5) * 0.5;
//         finalColor.rgb *= vignette;

//         // Enhance contrast slightly
//         finalColor.rgb = pow(finalColor.rgb, vec3(1.1));
//     }

//     if (finalColor.a == 0.0)
//         discard;

//     fragColor = clamp(finalColor, 0.0, 1.0);
// }


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

// Fluid style: 0=gray, 1=blue water
uniform int fluidStyle = 1;  // Default to blue water

// Color map types - extended with new gradient effects
uniform int colorMapType = 0; // 0: default, 1: blue, 2: purple, 3: cyan-yellow, 4: orange-grey
                             // 5: fire gradient, 6: ocean gradient, 7: plasma gradient
                             // 8: rainbow, 9: aurora, 10: lava

uniform int renderMode = 0; // 0: volume, 1: shell

// Define the position of the box - increased from 0.5 to 2.0 for larger display
const vec3 boxMin = vec3(-2.0);
const vec3 boxMax = vec3(2.0);

// Set the step parameters in ray marching
// If the simulation is not detailed enough, we should decrease stepSize or increase maxSteps
uniform int size;
float voxelLength = 1.0 / float(size);

// Simple hash function for noise generation
float hash(float n) {
    return fract(sin(n) * 43758.5453);
}

// 3D noise function for natural randomness
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

// Fractal Brownian Motion for more natural noise
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
    float stepSize = 0.4 / float(size);  // Reduced from 0.5 for more detailed sampling
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

        // Get cell-specific color type (if any)
        float cellColorType = texture(colorTex, texCoord).r;
        int localColorType = int(cellColorType + 0.5); // Round to nearest integer

        // Set color and alpha based on fluid style and color mapping
        vec3 color;
        float alphaValue;
        float depth = texCoord.y;  // Use y-coordinate for depth-based effects

        if (fluidStyle == 0) {  // Gray (original)
            color = vec3(1.0 - density);
            alphaValue = density * 0.15;
        } else if (fluidStyle == 1) {  // Blue water with color mapping
            // More precise color type determination
            int effectiveColorType;
            if (localColorType > 0) {
                // If cell has its own color, use it
                effectiveColorType = localColorType;
            } else if (density < 0.05) {
                // For very low density areas, use default coloring
                effectiveColorType = 0;
            } else {
                // For cells with no color but sufficient density, use global color
                effectiveColorType = colorMapType;
            }

            if (localColorType == 999) {  // Special value for obstacles
                // Use a fixed color for obstacles
                color = vec3(0.4, 0.4, 0.4);  // Simple gray color
                alphaValue = density * 0.8;   // More opaque
            }
            else if (effectiveColorType == 0) {
                // Original blue water color logic with enhanced depth effect
                float depthEffect = pow(depth, 1.2); // Makes depth effect more pronounced
                color = mix(waterDeepColor, waterBaseColor, depthEffect);

                // Add foam effect to high density areas
                color = mix(color, waterHighlightColor, pow(density, 2.5) * 0.6);

                // Add surface reflection at top of fluid with stronger effect
                float surfaceFactor = smoothstep(0.75, 1.0, texCoord.y);
                color = mix(color, waterSurfaceColor, surfaceFactor * 0.6);

                // Add edge highlights with more intensity
                color = mix(color, waterHighlightColor, clamp(edge * 1.2, 0.0, 0.6));

                // Adjust alpha value based on render mode
                if (renderMode == 1) {  // Shell rendering mode
                    alphaValue = density * 0.5 * (1.0 + 0.6 * depth);
                } else {
                    alphaValue = density * 0.18 * (1.0 + 0.5 * depth);
                }
            }
            else if (effectiveColorType == 1) {
                // Cyan-blue theme with enhanced contrast
                color = mix(vec3(0.0, 0.9, 1.0), vec3(0.0, 0.1, 0.6), pow(1.0-depth, 1.5));
                // Add highlights at density variations
                color = mix(color, vec3(0.6, 0.9, 1.0), pow(density, 3.0) * 0.7);
                alphaValue = density * 0.22;

                // Enhance shell mode rendering
                if (renderMode == 1) {
                    alphaValue *= 2.0;
                    color = mix(color, vec3(0.4, 0.8, 1.0), 0.2);
                }
            }
            else if (effectiveColorType == 2) {
                // Purple electric theme with more dramatic gradient
                color = mix(vec3(0.7, 0.0, 1.0), vec3(0.1, 0.0, 0.4), pow(1.0-depth, 1.3));
                // Add highlights
                color = mix(color, vec3(0.9, 0.5, 1.0), pow(density, 2.8) * 0.6);
                alphaValue = density * 0.22;

                // Enhance shell mode rendering
                if (renderMode == 1) {
                    alphaValue *= 2.0;
                    color = mix(color, vec3(0.8, 0.3, 1.0), 0.2);
                }
            }
            else if (effectiveColorType == 3) {
                // Cyan-yellow gradient with stronger contrast
                float gradientPos = depth * 0.7 + density * 0.3;
                color = mix(vec3(0.0, 0.9, 0.9), vec3(1.0, 0.9, 0.1), pow(gradientPos, 1.2));
                alphaValue = density * 0.22;

                // Enhance shell mode rendering
                if (renderMode == 1) {
                    alphaValue *= 2.0;
                    color = mix(color, vec3(0.5, 1.0, 0.5), 0.2);
                }
            }
            else if (effectiveColorType == 4) {
                // Orange-grey explosion effect with more contrast
                float gradientPos = depth * 0.6 + density * 0.4;
                color = mix(vec3(1.0, 0.4, 0.0), vec3(0.6, 0.6, 0.6), pow(gradientPos, 1.1));
                alphaValue = density * 0.22;

                // Enhance shell mode rendering
                if (renderMode == 1) {
                    alphaValue *= 2.0;
                    color = mix(color, vec3(1.0, 0.6, 0.2), 0.2);
                }
            }
            else if (effectiveColorType == 5) {
                // Fire gradient with natural irregular layers

                // Create a noise-based distortion of the depth
                // This creates irregular, flowing boundaries instead of straight lines
                float noiseScale = 5.0; // Controls the size of noise features
                float noiseStrength = 0.2; // Controls how much the noise distorts the layers

                // Generate base noise - varies in all three dimensions
                float baseNoise = fbm(vec3(texCoord.x * noiseScale,
                                          texCoord.y * noiseScale * 0.5,
                                          texCoord.z * noiseScale));

                // Create more detailed noise for fine details
                float detailNoise = noise(vec3(texCoord.x * noiseScale * 2.0,
                                              texCoord.y * noiseScale * 3.0,
                                              texCoord.z * noiseScale * 2.0));

                // Combine noise with depth to create irregular layering
                // Use both position and density to influence the layer position
                float layerPos = depth + (baseNoise * 0.3 + detailNoise * 0.1) * noiseStrength;
                // Add density influence - denser regions burn hotter (more yellow/white)
                layerPos += density * 0.15;

                // Create natural fire color gradient
                vec3 flameColor;
                if (layerPos < 0.25) {
                    // Inner core - bright white-yellow
                    flameColor = mix(vec3(1.0, 0.9, 0.5), vec3(1.0, 0.8, 0.0), smoothstep(0.0, 0.25, layerPos));
                } else if (layerPos < 0.5) {
                    // Mid flame - yellow to orange
                    flameColor = mix(vec3(1.0, 0.8, 0.0), vec3(1.0, 0.5, 0.0), smoothstep(0.25, 0.5, layerPos));
                } else if (layerPos < 0.75) {
                    // Outer flame - orange to red
                    flameColor = mix(vec3(1.0, 0.5, 0.0), vec3(0.9, 0.1, 0.0), smoothstep(0.5, 0.75, layerPos));
                } else {
                    // Fading edges - red to dark red
                    flameColor = mix(vec3(0.9, 0.1, 0.0), vec3(0.4, 0.0, 0.0), smoothstep(0.75, 1.0, layerPos));
                }

                // Add flickering effect
                float flicker = noise(vec3(texCoord.x, texCoord.z, time * 2.0)) * 0.2 + 0.9;
                flameColor *= flicker;

                // Add glow based on density and position
                float glowStrength = 0.4;
                vec3 glowColor = vec3(1.0, 0.3, 0.0) * density * glowStrength;
                flameColor += glowColor * (1.0 - layerPos); // More glow in the center/bottom

                // Final fire color
                color = flameColor;

                // Higher opacity for better visibility
                alphaValue = density * 0.28;

                // Enhance edge boundaries for more natural fire look
                if (edge > 0.1) {
                    // Add bright highlight at turbulent edges
                    color = mix(color, vec3(1.0, 0.8, 0.3), min(edge * 0.6, 0.5));
                }

                // Shell rendering enhancements
                if (renderMode == 1) {
                    alphaValue *= 2.0;
                    // Add highlighting in shell mode
                    color = mix(color, vec3(1.0, 0.7, 0.3), 0.15);
                }
            }
            else if (effectiveColorType == 6) {
                // ENHANCED Ocean gradient with stronger color contrast

                // Create noise-based distortion of the depth with increased strength
                float noiseScale = 4.0;
                float noiseStrength = 0.35; // Increased from 0.25 for more distortion

                // Generate base noise for large flowing patterns
                float baseNoise = fbm(vec3(texCoord.x * noiseScale,
                                          texCoord.y * noiseScale * 0.3,
                                          texCoord.z * noiseScale));

                // Create more detailed noise for small currents and eddies
                float detailNoise = noise(vec3(texCoord.x * noiseScale * 3.0,
                                              texCoord.y * noiseScale * 2.0,
                                              texCoord.z * noiseScale * 3.0));

                // Combine noise with depth for irregular layering
                // Create flowing water effect - stronger horizontal variation
                float layerPos = depth + (baseNoise * 0.5 + detailNoise * 0.2) * noiseStrength;
                // Add slight density influence - denser regions appear deeper
                layerPos -= density * 0.15;

                // Create more vibrant ocean color gradient with extreme contrast
                vec3 oceanColor;
                if (layerPos < 0.25) {
                    // Deep ocean - much darker blue to create stronger contrast
                    oceanColor = vec3(0.0, 0.0, 0.2);
                } else if (layerPos < 0.5) {
                    // Mid-depths - rich royal blue
                    oceanColor = mix(vec3(0.0, 0.0, 0.4), vec3(0.0, 0.3, 0.8), smoothstep(0.25, 0.5, layerPos));
                } else if (layerPos < 0.75) {
                    // Upper depths - vivid teal/turquoise with higher saturation
                    oceanColor = mix(vec3(0.0, 0.3, 0.8), vec3(0.0, 0.6, 1.0), smoothstep(0.5, 0.75, layerPos));
                } else {
                    // Surface waters - bright cyan to create more visible layers
                    oceanColor = mix(vec3(0.0, 0.6, 1.0), vec3(0.3, 0.8, 1.0), smoothstep(0.75, 1.0, layerPos));
                }

                // Add more dramatic underwater caustics effect
                float causticScale = 20.0;
                float caustics = sin(texCoord.x * causticScale + time * 1.2) *
                                sin(texCoord.z * causticScale + time * 1.8) * 0.5 + 0.5;
                // Stronger caustics effect
                oceanColor += vec3(0.0, 0.2, 0.3) * caustics * density * (layerPos + 0.3) * 1.5;

                // Add slight color variation based on density with stronger effect
                oceanColor *= (0.7 + density * 0.8);

                // Final ocean color
                color = oceanColor;

                // Higher opacity for better visibility
                alphaValue = density * 0.3;

                // Enhance edges for flowing water effect with stronger highlights
                if (edge > 0.1) {
                    // Add bright highlights at current edges
                    color = mix(color, vec3(0.7, 0.9, 1.0), min(edge * 0.7, 0.5));
                }

                // Shell rendering enhancements
                if (renderMode == 1) {
                    alphaValue *= 2.0;
                    // Add highlighting in shell mode
                    color = mix(color, vec3(0.4, 0.7, 1.0), 0.2);
                }
            }
            else if (effectiveColorType == 7) {
                // Plasma gradient with swirling electric patterns

                // Create noise-based distortion with higher frequency for electric look
                float noiseScale = 6.0;
                float noiseStrength = 0.3;

                // Generate base noise for large plasma patterns
                float baseNoise = fbm(vec3(texCoord.x * noiseScale,
                                         texCoord.y * noiseScale * 0.4,
                                         texCoord.z * noiseScale));

                // Create more detailed noise for electric arcs and tendrils
                float detailNoise = noise(vec3(texCoord.x * noiseScale * 4.0,
                                             texCoord.y * noiseScale * 3.0,
                                             texCoord.z * noiseScale * 4.0));

                // Combine noise with depth for swirling effect
                // Create electric plasma effect with more variation
                float layerPos = depth + (baseNoise * 0.35 + detailNoise * 0.2) * noiseStrength;
                // Add density influence - denser regions have more electric activity
                layerPos += density * 0.2;

                // Create electric plasma color gradient
                vec3 plasmaColor;
                if (layerPos < 0.25) {
                    // Core plasma - deep purple to violet
                    plasmaColor = mix(vec3(0.3, 0.0, 0.5), vec3(0.5, 0.0, 0.8), smoothstep(0.0, 0.25, layerPos));
                } else if (layerPos < 0.5) {
                    // Mid plasma - violet to magenta
                    plasmaColor = mix(vec3(0.5, 0.0, 0.8), vec3(0.9, 0.0, 0.9), smoothstep(0.25, 0.5, layerPos));
                } else if (layerPos < 0.75) {
                    // Outer plasma - magenta to hot orange
                    plasmaColor = mix(vec3(0.9, 0.0, 0.9), vec3(1.0, 0.3, 0.0), smoothstep(0.5, 0.75, layerPos));
                } else {
                    // Dissipating plasma - orange to yellow
                    plasmaColor = mix(vec3(1.0, 0.3, 0.0), vec3(1.0, 0.8, 0.0), smoothstep(0.75, 1.0, layerPos));
                }

                // Add electric arcing effect
                float arcNoise = noise(vec3(texCoord.x * 15.0, texCoord.y * 15.0, texCoord.z * 15.0 + time));
                float arcIntensity = pow(arcNoise, 5.0) * 2.0; // Makes arcs sharper and less frequent
                plasmaColor += vec3(0.8, 0.4, 1.0) * arcIntensity * density;

                // Add pulsing glow effect
                float pulseRate = 2.0; // Speed of pulsing
                float pulseNoise = noise(vec3(texCoord.x * 2.0, texCoord.y * 2.0, texCoord.z * 2.0 + time * pulseRate));
                plasmaColor *= 0.8 + pulseNoise * 0.4;

                // Final plasma color
                color = plasmaColor;

                alphaValue = density * 0.28;

                // Enhance edges for electric arc effect
                if (edge > 0.1) {
                    // Add bright highlight at energy boundaries
                    color = mix(color, vec3(0.9, 0.6, 1.0), min(edge * 0.7, 0.6));
                }

                // Shell rendering enhancements
                if (renderMode == 1) {
                    alphaValue *= 2.0;
                    // Add highlighting in shell mode
                    color = mix(color, vec3(0.8, 0.3, 1.0), 0.2);
                }
            }
            // NEW COLORMAP: Rainbow
            else if (effectiveColorType == 8) {
                // Rainbow gradient with flowing bands of color

                // Create noise distortion for irregular rainbow bands
                float noiseScale = 5.0;
                float noiseStrength = 0.25;

                // Base and detail noise
                float baseNoise = fbm(vec3(texCoord.x * noiseScale,
                                         texCoord.y * noiseScale * 0.3,
                                         texCoord.z * noiseScale));
                float detailNoise = noise(vec3(texCoord.x * noiseScale * 2.0,
                                            texCoord.z * noiseScale * 2.0,
                                            time * 0.2));

                // Combine with depth for flowing, irregular bands
                float bandPos = depth + (baseNoise * 0.3 + detailNoise * 0.1) * noiseStrength;
                // Add density influence
                bandPos += density * 0.1;

                // Prevent blue at the top by making sure bandPos never equals 1.0
                bandPos = min(bandPos, 0.99);

                // Create rainbow color spectrum
                vec3 rainbowColor;
                if (bandPos < 0.16) {
                    // Red
                    rainbowColor = vec3(1.0, 0.0, 0.0);
                } else if (bandPos < 0.33) {
                    // Orange
                    rainbowColor = mix(vec3(1.0, 0.0, 0.0), vec3(1.0, 0.5, 0.0), (bandPos - 0.16) / 0.17);
                } else if (bandPos < 0.5) {
                    // Yellow
                    rainbowColor = mix(vec3(1.0, 0.5, 0.0), vec3(1.0, 1.0, 0.0), (bandPos - 0.33) / 0.17);
                } else if (bandPos < 0.66) {
                    // Green
                    rainbowColor = mix(vec3(1.0, 1.0, 0.0), vec3(0.0, 0.8, 0.0), (bandPos - 0.5) / 0.16);
                } else if (bandPos < 0.83) {
                    // Blue
                    rainbowColor = mix(vec3(0.0, 0.8, 0.0), vec3(0.0, 0.0, 1.0), (bandPos - 0.66) / 0.17);
                } else {
                    // Violet
                    rainbowColor = mix(vec3(0.0, 0.0, 1.0), vec3(0.6, 0.0, 0.8), (bandPos - 0.83) / 0.17);
                }

                // Add shimmering effect
                float shimmer = noise(vec3(texCoord.x * 20.0, texCoord.z * 20.0, time * 3.0));
                rainbowColor += vec3(shimmer) * 0.1;

                // Adjust brightness based on density
                rainbowColor *= (0.7 + density * 0.5);

                // Final color
                color = rainbowColor;

                // Higher opacity
                alphaValue = density * 0.3;

                // Edge highlights
                if (edge > 0.1) {
                    color = mix(color, vec3(1.0), min(edge * 0.5, 0.4));
                }

                // Shell rendering enhancements
                if (renderMode == 1) {
                    alphaValue *= 2.0;
                    color = mix(color, vec3(1.0, 1.0, 1.0), 0.1);
                }
            }
            // NEW COLORMAP: Aurora
            else if (effectiveColorType == 9) {
                // Aurora borealis effect with dancing lights

                // Create time-varying noise for aurora movement
                float timeScale = time * 0.3;
                float noiseScale = 4.0;

                // Create flowing noise patterns
                float flowNoise = fbm(vec3(texCoord.x * noiseScale,
                                         texCoord.z * noiseScale,
                                         timeScale));
                float waveNoise = sin(texCoord.x * 10.0 + timeScale) *
                                 cos(texCoord.z * 8.0 + timeScale * 1.2) * 0.5 + 0.5;

                // Create vertical bands with distortion
                float bandPos = texCoord.y + (flowNoise * 0.3 + waveNoise * 0.1) * 0.4;
                // Density affects the band position
                bandPos += density * 0.1;

                // Auroral colors - greens, teals, purples
                vec3 auroraColor;
                if (bandPos < 0.3) {
                    // Deep purple to teal
                    auroraColor = mix(vec3(0.2, 0.0, 0.4), vec3(0.0, 0.4, 0.4), smoothstep(0.0, 0.3, bandPos));
                } else if (bandPos < 0.6) {
                    // Teal to bright green
                    auroraColor = mix(vec3(0.0, 0.4, 0.4), vec3(0.0, 1.0, 0.3), smoothstep(0.3, 0.6, bandPos));
                } else if (bandPos < 0.8) {
                    // Green to light cyan
                    auroraColor = mix(vec3(0.0, 1.0, 0.3), vec3(0.4, 1.0, 0.8), smoothstep(0.6, 0.8, bandPos));
                } else {
                    // Cyan to purple for high regions
                    auroraColor = mix(vec3(0.4, 1.0, 0.8), vec3(0.6, 0.3, 1.0), smoothstep(0.8, 1.0, bandPos));
                }

                // Add dancing light rays effect
                float rayIntensity = pow(sin(texCoord.x * 5.0 + timeScale * 2.0) *
                                       sin(texCoord.z * 5.0 + timeScale * 1.5), 2.0);
                auroraColor += vec3(0.2, 0.5, 0.2) * rayIntensity * (1.0 - depth) * 0.5;

                // Add subtle starfield effect for high regions
                if (depth > 0.7) {
                    float stars = step(0.98, noise(vec3(texCoord.x * 50.0, texCoord.y * 50.0, texCoord.z * 50.0)));
                    auroraColor += vec3(stars) * 0.5 * (depth - 0.7) / 0.3;
                }

                // Final color with brightness adjustment
                color = auroraColor * (0.8 + density * 0.5);

                // Higher opacity for better visibility
                alphaValue = density * 0.3;

                // Edge highlights
                if (edge > 0.15) {
                    color = mix(color, vec3(0.6, 1.0, 0.8), min(edge * 0.5, 0.4));
                }

                // Shell rendering enhancements
                if (renderMode == 1) {
                    alphaValue *= 2.0;
                    color = mix(color, vec3(0.0, 1.0, 0.5), 0.15);
                }
            }
            // NEW COLORMAP: Lava - Simplified version without simplex noise
            else if (effectiveColorType == 10) {
                // Lava flow with molten rock and cooling surface

                // Create noise for uneven lava flow
                float noiseScale = 3.0;
                float noiseStrength = 0.3;

                // Generate churning lava patterns
                float baseNoise = fbm(vec3(texCoord.x * noiseScale,
                                         texCoord.y * noiseScale * 0.2,
                                         texCoord.z * noiseScale));
                float detailNoise = noise(vec3(texCoord.x * noiseScale * 2.0,
                                            texCoord.y * noiseScale * 1.5,
                                            texCoord.z * noiseScale * 2.0 + time * 0.2));

                // Create bubbling, flowing lava
                float lavaPos = depth + (baseNoise * 0.4 + detailNoise * 0.2) * noiseStrength;
                // Denser regions are hotter
                lavaPos += density * 0.2;

                // Create molten lava color gradient
                vec3 lavaColor;
                if (lavaPos < 0.3) {
                    // Molten core - bright orange-yellow
                    lavaColor = mix(vec3(1.0, 0.8, 0.0), vec3(1.0, 0.5, 0.0), smoothstep(0.0, 0.3, lavaPos));
                } else if (lavaPos < 0.6) {
                    // Mid lava - orange to red
                    lavaColor = mix(vec3(1.0, 0.5, 0.0), vec3(0.9, 0.1, 0.0), smoothstep(0.3, 0.6, lavaPos));
                } else if (lavaPos < 0.85) {
                    // Cooling lava - red to dark red
                    lavaColor = mix(vec3(0.9, 0.1, 0.0), vec3(0.5, 0.0, 0.0), smoothstep(0.6, 0.85, lavaPos));
                } else {
                    // Solidified crust - dark red to black
                    lavaColor = mix(vec3(0.5, 0.0, 0.0), vec3(0.1, 0.0, 0.0), smoothstep(0.85, 1.0, lavaPos));
                }

                // Add bubbling/glowing effect using simplified approach
                float bubbleNoise = noise(vec3(texCoord.x * 15.0, texCoord.y * 10.0, texCoord.z * 15.0 + time * 0.5));
                float bubbleIntensity = bubbleNoise * bubbleNoise; // Squaring gives similar effect to step function
                lavaColor += vec3(1.0, 0.5, 0.0) * bubbleIntensity * density * 0.5;

                // Add simple cracks
                if (lavaPos > 0.7) {
                    float crackNoise = noise(vec3(texCoord.x * 20.0, texCoord.y * 5.0, texCoord.z * 20.0));
                    if (crackNoise > 0.7) {
                        lavaColor += vec3(1.0, 0.4, 0.0) * (lavaPos - 0.7) * 3.0;
                    }
                }

                // Final lava color with glow
                color = lavaColor;

                // Add emissive glow based on temperature (higher near bottom)
                color += vec3(0.5, 0.1, 0.0) * density * pow(1.0 - lavaPos, 2.0) * 0.3;

                // Higher opacity for better visibility
                alphaValue = density * 0.35;

                // Edge highlights for lava flow
                if (edge > 0.1) {
                    color = mix(color, vec3(1.0, 0.7, 0.0), min(edge * 0.6, 0.5));
                }

                // Shell rendering enhancements
                if (renderMode == 1) {
                    alphaValue *= 2.0;
                    color = mix(color, vec3(1.0, 0.3, 0.0), 0.2);
                }
            }
            else {
                // Fallback to default coloring
                color = waterBaseColor;
                alphaValue = density * 0.2;
            }

            // Enhanced edge highlights for all color schemes
            color = mix(color, vec3(1.0), clamp(edge * 1.5, 0.0, 0.7));

            // Special top surface handling to prevent blue artifacts
            if (texCoord.y > 0.95 && effectiveColorType > 0) {
                // Near the top surface, gradually fade to the scheme's own color instead of default blue
                float topFactor = (texCoord.y - 0.95) / 0.05;
                // Reduce any blue component to prevent blue artifacts
                color.b = mix(color.b, min(color.b, 0.5), topFactor);
            }

            // Add stronger Fresnel effect at glancing angles only for default water
            if (effectiveColorType == 0) {
                float fresnelEffect = fresnel(abs(dot(rayDir, vec3(0.0, 1.0, 0.0))), 0.05);
                color = mix(color, vec3(1.0), fresnelEffect * 0.4);
            }

            // Boost alpha in shell rendering mode
            if (renderMode == 1) {
                alphaValue *= 1.5;
            }
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
