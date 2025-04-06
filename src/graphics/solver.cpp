#include "solver.h"

int IX(int x, int y, int z, int N) {
    return x + N * (y + N * z);
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
    };
}
