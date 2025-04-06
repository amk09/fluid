#include "solver.h"

int IX(int x, int y, int z, int N) {
    return x + N * (y + N * z);
}

//For velocity fields, the component perpendicular to a wall is reflected.
//For scalar fields, values at the boundary are copied from adjacent cells.
void set_bnd(int N, int b, vector<float> &dataWrittenTo){
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
}


void diffuse_velocity(int b, vector<float> &dataWrittenTo, vector<float> &dataReadFrom, float diffuse, float dt){
    int size = dataWrittenTo.size();
    int N = std::cbrt(size);
    float a = dt * diffuse * size;

    for (int k = 0; k<20; k++){
        for(int x = 1; x<size-1; x++){
            for(int y=1; y<size-1; y++){
                for(int z = 1; z<size-2; z++){
                    int idx = IX(x, y, z, N);
                    dataWrittenTo[idx] =  (dataReadFrom[idx] + diffuse*(dataWrittenTo[IX(x-1, y, z, N)] + dataWrittenTo[IX(x+1, y, z, N)] +
                                                                        dataWrittenTo[IX(x, y-1, z, N)] + dataWrittenTo[IX(x, y+1, z, N)] +
                                                                        dataWrittenTo[IX(x, y, z-1, N)] + dataWrittenTo[IX(x, y, z+1, N)]
                                                                        )) / (1.f + 6.f*diffuse);
                }
            }
        }
        
    }
    set_bnd(N, b, dataWrittenTo);

}

void diffuse_density(int b, vector<float> &dataWrittenTo, vector<float> &dataReadFrom, float diffuse, float dt){
    int size = dataWrittenTo.size();
    int N = std::cbrt(size);
    float a = dt * diffuse * size;

    for (int k = 0; k<20; k++){
        for(int x = 1; x<size-1; x++){
            for(int y=1; y<size-1; y++){
                for(int z = 1; z<size-2; z++){
                    int idx = IX(x, y, z, N);
                    dataWrittenTo[idx] =  (dataReadFrom[idx] + diffuse*(dataWrittenTo[IX(x-1, y, z, N)] + dataWrittenTo[IX(x+1, y, z, N)] +
                                                                        dataWrittenTo[IX(x, y-1, z, N)] + dataWrittenTo[IX(x, y+1, z, N)] +
                                                                        dataWrittenTo[IX(x, y, z-1, N)] + dataWrittenTo[IX(x, y, z+1, N)]
                                                                        )) / (1.f + 6.f*diffuse);
                }
            }
        }
    }
    set_bnd(N, b, dataWrittenTo);
}

void advect(int b, vector<float> &dataWrittenTo, vector<float> &dataReadFrom, vector<float> &velocityX, vector<float> &velocityY, vector<float> &velocityZ, float dt){
    int size = dataWrittenTo.size();
    int N = std::round(std::cbrt(size));
    float dt0 = dt * N;
    for(int i = 1; i<N-1; i++){
        for(int j = 1; j < N-1; j++){
            for(int k = 1; k < N-1; k++){
                int idx = IX(i, j, k, N);

                //Backtrace particle position
                float x = i- dt0* velocityX[idx];
                float y = i- dt0* velocityY[idx];
                float z = i- dt0* velocityZ[idx];

                x = std::clamp(x, 0.5f, N - 1.5f);
                y = std::clamp(y, 0.5f, N - 1.5f);
                z = std::clamp(z, 0.5f, N - 1.5f);

                int i0 = (int)x;
                int j0 = (int)y;
                int k0 = (int)z;

                int i1 = i0 + 1;
                int j1 = j0 + 1;
                int k1 = k0 + 1;

                float s1 = x - i0, s0 = 1 - s1;
                float t1 = y - j0, t0 = 1 - t1;
                float u1 = z - k0, u0 = 1 - u1;

                dataWrittenTo[idx] =
                    s0 * (t0 * (u0 * dataReadFrom[IX(i0, j0, k0, N)] + u1 * dataReadFrom[IX(i0, j0, k1, N)]) +
                          t1 * (u0 * dataReadFrom[IX(i0, j1, k0, N)] + u1 * dataReadFrom[IX(i0, j1, k1, N)])) +
                    s1 * (t0 * (u0 * dataReadFrom[IX(i1, j0, k0, N)] + u1 * dataReadFrom[IX(i1, j0, k1, N)]) +
                          t1 * (u0 * dataReadFrom[IX(i1, j1, k0, N)] + u1 * dataReadFrom[IX(i1, j1, k1, N)]));


            }
        }
    }
    set_bnd(N, b, dataWrittenTo);
}

