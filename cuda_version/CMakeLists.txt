cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(fluidCuda CUDA CXX)

# Set CUDA specific flags
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --extended-lambda")

# Find required packages
find_package(CUDA REQUIRED)

# Add executable
add_executable(fluid_sim 
    src/main.cpp
    src/fluid_simulation.cu
)

# Set CUDA specific properties
set_target_properties(fluid_sim PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
)

# Include directories
target_include_directories(fluid_sim PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Link libraries
target_link_libraries(fluid_sim PRIVATE
    cuda
    cudart
) 