# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.25

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/asingir1/cs0320/fluidCuda

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/asingir1/cs0320/fluidCuda/build

# Include any dependencies generated for this target.
include CMakeFiles/fluid_sim.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/fluid_sim.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/fluid_sim.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/fluid_sim.dir/flags.make

CMakeFiles/fluid_sim.dir/src/main.cpp.o: CMakeFiles/fluid_sim.dir/flags.make
CMakeFiles/fluid_sim.dir/src/main.cpp.o: /home/asingir1/cs0320/fluidCuda/src/main.cpp
CMakeFiles/fluid_sim.dir/src/main.cpp.o: CMakeFiles/fluid_sim.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/asingir1/cs0320/fluidCuda/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/fluid_sim.dir/src/main.cpp.o"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/fluid_sim.dir/src/main.cpp.o -MF CMakeFiles/fluid_sim.dir/src/main.cpp.o.d -o CMakeFiles/fluid_sim.dir/src/main.cpp.o -c /home/asingir1/cs0320/fluidCuda/src/main.cpp

CMakeFiles/fluid_sim.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/fluid_sim.dir/src/main.cpp.i"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/asingir1/cs0320/fluidCuda/src/main.cpp > CMakeFiles/fluid_sim.dir/src/main.cpp.i

CMakeFiles/fluid_sim.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/fluid_sim.dir/src/main.cpp.s"
	/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/asingir1/cs0320/fluidCuda/src/main.cpp -o CMakeFiles/fluid_sim.dir/src/main.cpp.s

CMakeFiles/fluid_sim.dir/src/fluid_simulation.cu.o: CMakeFiles/fluid_sim.dir/flags.make
CMakeFiles/fluid_sim.dir/src/fluid_simulation.cu.o: CMakeFiles/fluid_sim.dir/includes_CUDA.rsp
CMakeFiles/fluid_sim.dir/src/fluid_simulation.cu.o: /home/asingir1/cs0320/fluidCuda/src/fluid_simulation.cu
CMakeFiles/fluid_sim.dir/src/fluid_simulation.cu.o: CMakeFiles/fluid_sim.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/asingir1/cs0320/fluidCuda/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CUDA object CMakeFiles/fluid_sim.dir/src/fluid_simulation.cu.o"
	/local/projects/cuda12/cuda12.2.2/bin/nvcc -forward-unknown-to-host-compiler $(CUDA_DEFINES) $(CUDA_INCLUDES) $(CUDA_FLAGS) -MD -MT CMakeFiles/fluid_sim.dir/src/fluid_simulation.cu.o -MF CMakeFiles/fluid_sim.dir/src/fluid_simulation.cu.o.d -x cu -rdc=true -c /home/asingir1/cs0320/fluidCuda/src/fluid_simulation.cu -o CMakeFiles/fluid_sim.dir/src/fluid_simulation.cu.o

CMakeFiles/fluid_sim.dir/src/fluid_simulation.cu.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CUDA source to CMakeFiles/fluid_sim.dir/src/fluid_simulation.cu.i"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_PREPROCESSED_SOURCE

CMakeFiles/fluid_sim.dir/src/fluid_simulation.cu.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CUDA source to assembly CMakeFiles/fluid_sim.dir/src/fluid_simulation.cu.s"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_ASSEMBLY_SOURCE

# Object files for target fluid_sim
fluid_sim_OBJECTS = \
"CMakeFiles/fluid_sim.dir/src/main.cpp.o" \
"CMakeFiles/fluid_sim.dir/src/fluid_simulation.cu.o"

# External object files for target fluid_sim
fluid_sim_EXTERNAL_OBJECTS =

CMakeFiles/fluid_sim.dir/cmake_device_link.o: CMakeFiles/fluid_sim.dir/src/main.cpp.o
CMakeFiles/fluid_sim.dir/cmake_device_link.o: CMakeFiles/fluid_sim.dir/src/fluid_simulation.cu.o
CMakeFiles/fluid_sim.dir/cmake_device_link.o: CMakeFiles/fluid_sim.dir/build.make
CMakeFiles/fluid_sim.dir/cmake_device_link.o: CMakeFiles/fluid_sim.dir/deviceLinkLibs.rsp
CMakeFiles/fluid_sim.dir/cmake_device_link.o: CMakeFiles/fluid_sim.dir/deviceObjects1
CMakeFiles/fluid_sim.dir/cmake_device_link.o: CMakeFiles/fluid_sim.dir/dlink.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/asingir1/cs0320/fluidCuda/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CUDA device code CMakeFiles/fluid_sim.dir/cmake_device_link.o"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fluid_sim.dir/dlink.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/fluid_sim.dir/build: CMakeFiles/fluid_sim.dir/cmake_device_link.o
.PHONY : CMakeFiles/fluid_sim.dir/build

# Object files for target fluid_sim
fluid_sim_OBJECTS = \
"CMakeFiles/fluid_sim.dir/src/main.cpp.o" \
"CMakeFiles/fluid_sim.dir/src/fluid_simulation.cu.o"

# External object files for target fluid_sim
fluid_sim_EXTERNAL_OBJECTS =

fluid_sim: CMakeFiles/fluid_sim.dir/src/main.cpp.o
fluid_sim: CMakeFiles/fluid_sim.dir/src/fluid_simulation.cu.o
fluid_sim: CMakeFiles/fluid_sim.dir/build.make
fluid_sim: CMakeFiles/fluid_sim.dir/cmake_device_link.o
fluid_sim: CMakeFiles/fluid_sim.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/asingir1/cs0320/fluidCuda/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable fluid_sim"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fluid_sim.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/fluid_sim.dir/build: fluid_sim
.PHONY : CMakeFiles/fluid_sim.dir/build

CMakeFiles/fluid_sim.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/fluid_sim.dir/cmake_clean.cmake
.PHONY : CMakeFiles/fluid_sim.dir/clean

CMakeFiles/fluid_sim.dir/depend:
	cd /home/asingir1/cs0320/fluidCuda/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/asingir1/cs0320/fluidCuda /home/asingir1/cs0320/fluidCuda /home/asingir1/cs0320/fluidCuda/build /home/asingir1/cs0320/fluidCuda/build /home/asingir1/cs0320/fluidCuda/build/CMakeFiles/fluid_sim.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/fluid_sim.dir/depend

