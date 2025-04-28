#include "solver.h"
#include "graphics/obstacle.h"
#include <omp.h>

int IX(int x, int y, int z, int N) {
    return x + N * (y + N * z);
}

void diffuse_velocity(int b, vector<float> &dataWrittenTo, vector<float> &dataReadFrom, float viscosity, float dt, int iter, int N){
    float a = dt * viscosity * (N - 2) * (N - 2);
    lin_solve(b, dataWrittenTo, dataReadFrom, a, 1.f + 6.f*viscosity, iter, N);
}

void diffuse_density(int b, vector<float> &dataWrittenTo, vector<float> &dataReadFrom, float diffuse, float dt, int iter, int N){
    float a = dt * diffuse * (N - 2) * (N - 2);
    lin_solve(b, dataWrittenTo, dataReadFrom, a, 1.f + 6.f*diffuse, iter, N);
}

void lin_solve(int b, vector<float> &dataWrittenTo, vector<float> &dataReadFrom, float a, float c, int iter, int N){
    float cRecip = 1.0f / c;

    for (int k = 0; k < iter; ++k){
        #pragma omp parallel for collapse(3)
        for (int z = 1; z < N - 1; ++z){
            for (int y = 1; y < N - 1; ++y){
                for (int x = 1; x < N - 1; ++x){
                    dataWrittenTo[IX(x, y, z, N)] =
                        (dataReadFrom[IX(x, y, z, N)] +
                         a * (dataWrittenTo[IX(x-1, y, z, N)] +
                              dataWrittenTo[IX(x+1, y, z, N)] +
                              dataWrittenTo[IX(x, y-1, z, N)] +
                              dataWrittenTo[IX(x, y+1, z, N)] +
                              dataWrittenTo[IX(x, y, z-1, N)] +
                              dataWrittenTo[IX(x, y, z+1, N)])) * cRecip;
                }
            }
        }
        set_bnd(b, dataWrittenTo, N);
    }
}

void advect(int b, vector<float> &dataWrittenTo, vector<float> &dataReadFrom, vector<float> &velocityX, vector<float> &velocityY, vector<float> &velocityZ, float dt, int N){
    
    float dt0 = dt * (N - 2);
    for (int k = 1; k < N - 1; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {

                int idx = IX(i, j, k, N);

                // Backtrace particle position
                float x = i - dt0 * velocityX[idx];
                float y = j - dt0 * velocityY[idx];
                float z = k - dt0 * velocityZ[idx];

                // The real particle position is from 0.5 to N - 0.5f
                x = std::clamp(x, 0.5f, N - 1.5f);
                y = std::clamp(y, 0.5f, N - 1.5f);
                z = std::clamp(z, 0.5f, N - 1.5f);

                int i0 = (int)x;
                int j0 = (int)y;
                int k0 = (int)z;

                int i1 = i0 + 1;
                int j1 = j0 + 1;
                int k1 = k0 + 1;

                float s1 = x - i0, s0 = 1.0f - s1;
                float t1 = y - j0, t0 = 1.0f - t1;
                float u1 = z - k0, u0 = 1.0f - u1;

                dataWrittenTo[idx] =
                    s0 * (t0 * (u0 * dataReadFrom[IX(i0, j0, k0, N)] + u1 * dataReadFrom[IX(i0, j0, k1, N)]) +
                          t1 * (u0 * dataReadFrom[IX(i0, j1, k0, N)] + u1 * dataReadFrom[IX(i0, j1, k1, N)])) +
                    s1 * (t0 * (u0 * dataReadFrom[IX(i1, j0, k0, N)] + u1 * dataReadFrom[IX(i1, j0, k1, N)]) +
                          t1 * (u0 * dataReadFrom[IX(i1, j1, k0, N)] + u1 * dataReadFrom[IX(i1, j1, k1, N)]));
            }
        }
    }

    set_bnd(b, dataWrittenTo, N);
}

void project(vector<float> &velocityX, vector<float> &velocityY, vector<float> &velocityZ, vector<float> &p, vector<float> div, int iter, int N){
    // Step 1: Compute divergence
    for (int k = 1; k < N - 1; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                int idx = IX(i, j, k, N);
                div[idx] = -0.5f * (
                               velocityX[IX(i + 1, j, k, N)] - velocityX[IX(i - 1, j, k, N)] +
                               velocityY[IX(i, j + 1, k, N)] - velocityY[IX(i, j - 1, k, N)] +
                               velocityZ[IX(i, j, k + 1, N)] - velocityZ[IX(i, j, k - 1, N)]
                               ) / N;

                p[idx] = 0.0f;
            }
        }
    }

    set_bnd(0, div, N);
    set_bnd(0, p, N);

    // Step 2: Solve pressure Poisson equation ∇²p = div
    lin_solve(0, p, div, 1, 6, iter, N);  // a=1, c=6 as per standard

    // Step 3: Subtract pressure gradient to get divergence-free velocity
    for (int k = 1; k < N - 1; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                int idx = IX(i, j, k, N);

                velocityX[idx] -= 0.5f * (p[IX(i + 1, j, k, N)] - p[IX(i - 1, j, k, N)]) * N;
                velocityY[idx] -= 0.5f * (p[IX(i, j + 1, k, N)] - p[IX(i, j - 1, k, N)]) * N;
                velocityZ[idx] -= 0.5f * (p[IX(i, j, k + 1, N)] - p[IX(i, j, k - 1, N)]) * N;
            }
        }
    }

    set_bnd(1, velocityX, N);
    set_bnd(2, velocityY, N);
    set_bnd(3, velocityZ, N);
}

//For velocity fields, the component perpendicular to a wall is reflected.
//For scalar fields, values at the boundary are copied from adjacent cells.
void set_bnd(int b, vector<float> &dataWrittenTo, int N){

    // New
    if (!g_obstacle.empty()) {
        for (int k = 1; k < N-1; k++) {
            for (int j = 1; j < N-1; j++) {
                for (int i = 1; i < N-1; i++) {
                    int idx = IX(i, j, k, N);
                    if (g_obstacle[idx] == 1) {
                        if (g_obstacle[IX(i-1,j,k,N)] == 0)
                            dataWrittenTo[IX(i-1,j,k,N)] = (b == 1) ? -dataWrittenTo[idx] : dataWrittenTo[idx];
                        if (g_obstacle[IX(i+1,j,k,N)] == 0)
                            dataWrittenTo[IX(i+1,j,k,N)] = (b == 1) ? -dataWrittenTo[idx] : dataWrittenTo[idx];
                        if (g_obstacle[IX(i,j-1,k,N)] == 0)
                            dataWrittenTo[IX(i,j-1,k,N)] = (b == 2) ? -dataWrittenTo[idx] : dataWrittenTo[idx];
                        if (g_obstacle[IX(i,j+1,k,N)] == 0)
                            dataWrittenTo[IX(i,j+1,k,N)] = (b == 2) ? -dataWrittenTo[idx] : dataWrittenTo[idx];
                        if (g_obstacle[IX(i,j,k-1,N)] == 0)
                            dataWrittenTo[IX(i,j,k-1,N)] = (b == 3) ? -dataWrittenTo[idx] : dataWrittenTo[idx];
                        if (g_obstacle[IX(i,j,k+1,N)] == 0)
                            dataWrittenTo[IX(i,j,k+1,N)] = (b == 3) ? -dataWrittenTo[idx] : dataWrittenTo[idx];
                    }
                }
            }
        }
    }

    for (int i = 1; i < N - 1; i++) {
        for (int j = 1; j < N - 1; j++) {
            //for Z axis
            dataWrittenTo[IX(i, j, 0, N)]     = (b == 3 ? -dataWrittenTo[IX(i, j, 1, N)]     : dataWrittenTo[IX(i, j, 1, N)]);
            dataWrittenTo[IX(i, j, N-1, N)]   = (b == 3 ? -dataWrittenTo[IX(i, j, N-2, N)]   : dataWrittenTo[IX(i, j, N-2, N)]);

            // For Y axis
            dataWrittenTo[IX(i, 0, j, N)]     = (b == 2 ? -dataWrittenTo[IX(i, 1, j, N)]     : dataWrittenTo[IX(i, 1, j, N)]);
            dataWrittenTo[IX(i, N-1, j, N)]   = (b == 2 ? -dataWrittenTo[IX(i, N-2, j, N)]   : dataWrittenTo[IX(i, N-2, j, N)]);

            // For X axis
            dataWrittenTo[IX(0, i, j, N)]     = (b == 1 ? -dataWrittenTo[IX(1, i, j, N)]     : dataWrittenTo[IX(1, i, j, N)]);
            dataWrittenTo[IX(N-1, i, j, N)]   = (b == 1 ? -dataWrittenTo[IX(N-2, i, j, N)]   : dataWrittenTo[IX(N-2, i, j, N)]);
        }
    }

    //optional: average neighboring inner cells
    dataWrittenTo[IX(0,     0,     0,     N)] = (dataWrittenTo[IX(1, 0, 0, N)] + dataWrittenTo[IX(0, 1, 0, N)] + dataWrittenTo[IX(0, 0, 1, N)]) / 3.0f;
    dataWrittenTo[IX(N - 1, 0,     0,     N)] = (dataWrittenTo[IX(N - 2, 0, 0, N)] + dataWrittenTo[IX(N - 1, 1, 0, N)] + dataWrittenTo[IX(N - 1, 0, 1, N)]) / 3.0f;
    dataWrittenTo[IX(0,     N - 1, 0,     N)] = (dataWrittenTo[IX(1, N - 1, 0, N)] + dataWrittenTo[IX(0, N - 2, 0, N)] + dataWrittenTo[IX(0, N - 1, 1, N)]) / 3.0f;
    dataWrittenTo[IX(0,     0,     N - 1, N)] = (dataWrittenTo[IX(1, 0, N - 1, N)] + dataWrittenTo[IX(0, 1, N - 1, N)] + dataWrittenTo[IX(0, 0, N - 2, N)]) / 3.0f;
    dataWrittenTo[IX(N - 1, N - 1, 0,     N)] = (dataWrittenTo[IX(N - 2, N - 1, 0, N)] + dataWrittenTo[IX(N - 1, N - 2, 0, N)] + dataWrittenTo[IX(N - 1, N - 1, 1, N)]) / 3.0f;
    dataWrittenTo[IX(N - 1, 0,     N - 1, N)] = (dataWrittenTo[IX(N - 2, 0, N - 1, N)] + dataWrittenTo[IX(N - 1, 1, N - 1, N)] + dataWrittenTo[IX(N - 1, 0, N - 2, N)]) / 3.0f;
    dataWrittenTo[IX(0,     N - 1, N - 1, N)] = (dataWrittenTo[IX(1, N - 1, N - 1, N)] + dataWrittenTo[IX(0, N - 2, N - 1, N)] + dataWrittenTo[IX(0, N - 1, N - 2, N)]) / 3.0f;
    dataWrittenTo[IX(N - 1, N - 1, N - 1, N)] = (dataWrittenTo[IX(N - 2, N - 1, N - 1, N)] + dataWrittenTo[IX(N - 1, N - 2, N - 1, N)] + dataWrittenTo[IX(N - 1, N - 1, N - 2, N)]) / 3.0f;

}


void conserveDensity(std::vector<float> &density, float targetDensity, int N) {
    float currentDensity = 0.0f;
    
    // Calculate current total density
    for (int k = 0; k < N; k++) {
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < N; i++) {
                currentDensity += density[IX(i, j, k, N)];
            }
        }
    }

    if (currentDensity < 0.000001f) {
        return;
    }
    
    float factor = targetDensity / currentDensity;
    
    for (int k = 0; k < N; k++) {
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < N; i++) {
                density[IX(i, j, k, N)] *= factor;
            }
        }
    }
}





void addFountainForce(std::vector<float> &velocityX, std::vector<float> &velocityY,
                      std::vector<float> &velocityZ, std::vector<float> &density,
                      int N, float strength, float dt) {
    int center = N / 2;
    int fountainRadius = N / 20;  // Reduced radius for more concentrated fountain
    float fountainStrength = strength * dt * 1.5f;  // Increased strength

    // Apply upward force in a focused area at bottom
    for (int k = center - fountainRadius; k <= center + fountainRadius; k++) {
        for (int i = center - fountainRadius; i <= center + fountainRadius; i++) {
            // Calculate distance from center
            float dx = i - center;
            float dz = k - center;
            float dist2 = dx*dx + dz*dz;

            if (dist2 > fountainRadius*fountainRadius)
                continue;

            // Apply force stronger at center, weaker at edges
            float radiusFactor = 1.0f - sqrt(dist2) / fountainRadius;
            float force = fountainStrength * radiusFactor * radiusFactor;

            // Apply to bottom cells
            for (int j = 0; j < 3; j++) {
                int idx = IX(i, j, k, N);

                // Add upward velocity
                velocityY[idx] += force * 1.5f;  // More upward force

                // Add very slight X/Z variation to avoid perfect symmetry
                velocityX[idx] += (((float)rand()/RAND_MAX) * 0.1f - 0.05f) * force * 0.1f;
                velocityZ[idx] += (((float)rand()/RAND_MAX) * 0.1f - 0.05f) * force * 0.1f;
            }
        }
    }

    set_bnd(1, velocityX, N);
    set_bnd(2, velocityY, N);
    set_bnd(3, velocityZ, N);
}


void circulateDensity(std::vector<float> &density, int N, float dt, float rate) {
    int center = N / 2;
    int sourceRadius = N / 20;  // Smaller source radius
    int sinkRadius = N / 15;    // Smaller sink radius
    float transferRate = rate * dt * 1.5f;  // Increased transfer rate

    std::vector<float> densityChange(density.size(), 0.0f);
    float totalRemoved = 0.0f;

    // First pass: remove density from top of fluid in a concentrated area
    for (int k = center - sinkRadius; k <= center + sinkRadius; k++) {
        for (int i = center - sinkRadius; i <= center + sinkRadius; i++) {
            float dx = i - center;
            float dz = k - center;
            float dist2 = dx*dx + dz*dz;

            if (dist2 > sinkRadius*sinkRadius)
                continue;

            // Remove from top cells
            for (int j = N-5; j < N-1; j++) {
                int idx = IX(i, j, k, N);
                float amountToRemove = density[idx] * transferRate;
                densityChange[idx] -= amountToRemove;
                totalRemoved += amountToRemove;
            }
        }
    }

    // Second pass: add density to the fountain source (concentrated)
    int fountainPoints = 0;
    for (int k = center - sourceRadius; k <= center + sourceRadius; k++) {
        for (int i = center - sourceRadius; i <= center + sourceRadius; i++) {
            float dx = i - center;
            float dz = k - center;
            float dist2 = dx*dx + dz*dz;

            if (dist2 > sourceRadius*sourceRadius)
                continue;

            // Count points for bottom source
            for (int j = 0; j < 3; j++) {
                fountainPoints++;
            }
        }
    }

    float addPerCell = fountainPoints > 0 ? totalRemoved / fountainPoints : 0;

    for (int k = center - sourceRadius; k <= center + sourceRadius; k++) {
        for (int i = center - sourceRadius; i <= center + sourceRadius; i++) {
            float dx = i - center;
            float dz = k - center;
            float dist2 = dx*dx + dz*dz;

            if (dist2 > sourceRadius*sourceRadius)
                continue;

            // Add to bottom cells
            for (int j = 0; j < 3; j++) {
                int idx = IX(i, j, k, N);
                densityChange[idx] += addPerCell;
            }
        }
    }

    // Apply density changes
    for (int i = 0; i < density.size(); i++) {
        density[i] += densityChange[i];
        // Clamp density to [0, 1]
        density[i] = std::max(0.0f, std::min(density[i], 1.0f));
    }
}

void addVorticityConfinement(std::vector<float> &vX, std::vector<float> &vY,
                             std::vector<float> &vZ, float dt, int N, float vorticityStrength) {
    // Temporary arrays for curl (vorticity) calculation
    std::vector<float> curlX(N*N*N, 0.0f);
    std::vector<float> curlY(N*N*N, 0.0f);
    std::vector<float> curlZ(N*N*N, 0.0f);
    std::vector<float> curlMagnitude(N*N*N, 0.0f);

// Step 1: Calculate curl (∇ × v) at each cell
#pragma omp parallel for collapse(3)
    for (int k = 1; k < N-1; k++) {
        for (int j = 1; j < N-1; j++) {
            for (int i = 1; i < N-1; i++) {
                int idx = IX(i, j, k, N);

                // Calculate partial derivatives using central differences
                float dvz_dy = (vZ[IX(i, j+1, k, N)] - vZ[IX(i, j-1, k, N)]) * 0.5f;
                float dvy_dz = (vY[IX(i, j, k+1, N)] - vY[IX(i, j, k-1, N)]) * 0.5f;

                float dvx_dz = (vX[IX(i, j, k+1, N)] - vX[IX(i, j, k-1, N)]) * 0.5f;
                float dvz_dx = (vZ[IX(i+1, j, k, N)] - vZ[IX(i-1, j, k, N)]) * 0.5f;

                float dvy_dx = (vY[IX(i+1, j, k, N)] - vY[IX(i-1, j, k, N)]) * 0.5f;
                float dvx_dy = (vX[IX(i, j+1, k, N)] - vX[IX(i, j-1, k, N)]) * 0.5f;

                // Curl components (vorticity)
                curlX[idx] = dvz_dy - dvy_dz;
                curlY[idx] = dvx_dz - dvz_dx;
                curlZ[idx] = dvy_dx - dvx_dy;

                // Calculate curl magnitude
                curlMagnitude[idx] = sqrtf(
                    curlX[idx] * curlX[idx] +
                    curlY[idx] * curlY[idx] +
                    curlZ[idx] * curlZ[idx]
                    );
            }
        }
    }

// Step 2: Calculate normalized vorticity confinement force: N × ω
// N is the normalized gradient of the curl magnitude
#pragma omp parallel for collapse(3)
    for (int k = 1; k < N-1; k++) {
        for (int j = 1; j < N-1; j++) {
            for (int i = 1; i < N-1; i++) {
                int idx = IX(i, j, k, N);

                // Skip cells with very small curl
                if (curlMagnitude[idx] < 0.000001f) continue;

                // Calculate gradient of curl magnitude using central differences
                float dx = (curlMagnitude[IX(i+1, j, k, N)] - curlMagnitude[IX(i-1, j, k, N)]) * 0.5f;
                float dy = (curlMagnitude[IX(i, j+1, k, N)] - curlMagnitude[IX(i, j-1, k, N)]) * 0.5f;
                float dz = (curlMagnitude[IX(i, j, k+1, N)] - curlMagnitude[IX(i, j, k-1, N)]) * 0.5f;

                // Normalize gradient
                float length = sqrtf(dx*dx + dy*dy + dz*dz);
                if (length > 0.000001f) {
                    dx /= length;
                    dy /= length;
                    dz /= length;

                    // Calculate force as cross product: N × ω (normalized gradient × curl)
                    float forceX = (dy * curlZ[idx] - dz * curlY[idx]) * vorticityStrength;
                    float forceY = (dz * curlX[idx] - dx * curlZ[idx]) * vorticityStrength;
                    float forceZ = (dx * curlY[idx] - dy * curlX[idx]) * vorticityStrength;

                    // Apply force scaled by time step
                    vX[idx] += forceX * dt;
                    vY[idx] += forceY * dt;
                    vZ[idx] += forceZ * dt;
                }
            }
        }
    }

    // Apply boundary conditions
    set_bnd(1, vX, N);
    set_bnd(2, vY, N);
    set_bnd(3, vZ, N);
}
