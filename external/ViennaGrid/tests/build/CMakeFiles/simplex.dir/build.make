# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Produce verbose output by default.
VERBOSE = 1

# Suppress display of executed commands.
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
RM = /usr/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/weinbub/git/ViennaGrid/tests

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/weinbub/git/ViennaGrid/tests/build

# Include any dependencies generated for this target.
include CMakeFiles/simplex.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/simplex.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/simplex.dir/flags.make

CMakeFiles/simplex.dir/src/simplex.cpp.o: CMakeFiles/simplex.dir/flags.make
CMakeFiles/simplex.dir/src/simplex.cpp.o: ../src/simplex.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/weinbub/git/ViennaGrid/tests/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/simplex.dir/src/simplex.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/simplex.dir/src/simplex.cpp.o -c /home/weinbub/git/ViennaGrid/tests/src/simplex.cpp

CMakeFiles/simplex.dir/src/simplex.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simplex.dir/src/simplex.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/weinbub/git/ViennaGrid/tests/src/simplex.cpp > CMakeFiles/simplex.dir/src/simplex.cpp.i

CMakeFiles/simplex.dir/src/simplex.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simplex.dir/src/simplex.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/weinbub/git/ViennaGrid/tests/src/simplex.cpp -o CMakeFiles/simplex.dir/src/simplex.cpp.s

CMakeFiles/simplex.dir/src/simplex.cpp.o.requires:
.PHONY : CMakeFiles/simplex.dir/src/simplex.cpp.o.requires

CMakeFiles/simplex.dir/src/simplex.cpp.o.provides: CMakeFiles/simplex.dir/src/simplex.cpp.o.requires
	$(MAKE) -f CMakeFiles/simplex.dir/build.make CMakeFiles/simplex.dir/src/simplex.cpp.o.provides.build
.PHONY : CMakeFiles/simplex.dir/src/simplex.cpp.o.provides

CMakeFiles/simplex.dir/src/simplex.cpp.o.provides.build: CMakeFiles/simplex.dir/src/simplex.cpp.o

# Object files for target simplex
simplex_OBJECTS = \
"CMakeFiles/simplex.dir/src/simplex.cpp.o"

# External object files for target simplex
simplex_EXTERNAL_OBJECTS =

simplex: CMakeFiles/simplex.dir/src/simplex.cpp.o
simplex: CMakeFiles/simplex.dir/build.make
simplex: CMakeFiles/simplex.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable simplex"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/simplex.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/simplex.dir/build: simplex
.PHONY : CMakeFiles/simplex.dir/build

CMakeFiles/simplex.dir/requires: CMakeFiles/simplex.dir/src/simplex.cpp.o.requires
.PHONY : CMakeFiles/simplex.dir/requires

CMakeFiles/simplex.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/simplex.dir/cmake_clean.cmake
.PHONY : CMakeFiles/simplex.dir/clean

CMakeFiles/simplex.dir/depend:
	cd /home/weinbub/git/ViennaGrid/tests/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/weinbub/git/ViennaGrid/tests /home/weinbub/git/ViennaGrid/tests /home/weinbub/git/ViennaGrid/tests/build /home/weinbub/git/ViennaGrid/tests/build /home/weinbub/git/ViennaGrid/tests/build/CMakeFiles/simplex.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/simplex.dir/depend

