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

# Utility rule file for ExperimentalConfigure.

CMakeFiles/ExperimentalConfigure:
	/usr/bin/ctest -D ExperimentalConfigure

ExperimentalConfigure: CMakeFiles/ExperimentalConfigure
ExperimentalConfigure: CMakeFiles/ExperimentalConfigure.dir/build.make
.PHONY : ExperimentalConfigure

# Rule to build all files generated by this target.
CMakeFiles/ExperimentalConfigure.dir/build: ExperimentalConfigure
.PHONY : CMakeFiles/ExperimentalConfigure.dir/build

CMakeFiles/ExperimentalConfigure.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ExperimentalConfigure.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ExperimentalConfigure.dir/clean

CMakeFiles/ExperimentalConfigure.dir/depend:
	cd /home/weinbub/git/ViennaGrid/tests/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/weinbub/git/ViennaGrid/tests /home/weinbub/git/ViennaGrid/tests /home/weinbub/git/ViennaGrid/tests/build /home/weinbub/git/ViennaGrid/tests/build /home/weinbub/git/ViennaGrid/tests/build/CMakeFiles/ExperimentalConfigure.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ExperimentalConfigure.dir/depend

