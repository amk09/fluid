#!/bin/bash

# # # Define Conda environment name
# CONDA_ENV_PATH="/users/jhan192/my_conda_env"

# # # Activate Conda environment
# echo "Activating Conda environment..."
# source ~/miniconda3/etc/profile.d/conda.sh  # Ensure Conda is loaded
# conda activate /users/jhan192/my_conda_env

# Set the number of threads for OpenMP
export OMP_NUM_THREADS=$(nproc)
export LD_LIBRARY_PATH=$CONDA_ENV_PATH/lib:$LD_LIBRARY_PATH

# Remove old build files
echo "🚀 Cleaning previous build..."
rm -rf build
mkdir build && cd build

# Run CMake to configure the project
echo "⚙️  Configuring CMake..."
cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS="-fopenmp" \
         -DCMAKE_EXE_LINKER_FLAGS="-L/users/jhan192/my_conda_env/lib -Wl,-rpath,/users/jhan192/my_conda_env/lib -lstdc++"

# Check if CMake succeeded
if [ $? -ne 0 ]; then
    echo "❌ CMake configuration failed!"
    exit 1
fi

# Compile the project using all available cores
echo "🔨 Building the project... with the number of $(nproc) cores"
make -j$(nproc)

# Check if the build succeeded
if [ $? -ne 0 ]; then
    echo "Build failed!"
    exit 1
fi
cd ../
# if [ ! -f path ]; then
#     echo "Error: Executable 'path' not found in build directory!"
#     exit 1
# fi

# Run the executable with the provided scene file
SCENE_FILE="template_inis/final/refraction.ini"
echo "Running the program with scene: $SCENE_FILE"
./build/path $SCENE_FILE
