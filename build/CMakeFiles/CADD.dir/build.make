# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.0

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

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

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/srinath/CADD_all/CADD_ben2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/srinath/CADD_all/CADD_ben2/build

# Include any dependencies generated for this target.
include CMakeFiles/CADD.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/CADD.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/CADD.dir/flags.make

CMakeFiles/CADD.dir/field.f.o: CMakeFiles/CADD.dir/flags.make
CMakeFiles/CADD.dir/field.f.o: ../field.f
	$(CMAKE_COMMAND) -E cmake_progress_report /home/srinath/CADD_all/CADD_ben2/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building Fortran object CMakeFiles/CADD.dir/field.f.o"
	/opt/intel/bin/ifort  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /home/srinath/CADD_all/CADD_ben2/field.f -o CMakeFiles/CADD.dir/field.f.o

CMakeFiles/CADD.dir/field.f.o.requires:
.PHONY : CMakeFiles/CADD.dir/field.f.o.requires

CMakeFiles/CADD.dir/field.f.o.provides: CMakeFiles/CADD.dir/field.f.o.requires
	$(MAKE) -f CMakeFiles/CADD.dir/build.make CMakeFiles/CADD.dir/field.f.o.provides.build
.PHONY : CMakeFiles/CADD.dir/field.f.o.provides

CMakeFiles/CADD.dir/field.f.o.provides.build: CMakeFiles/CADD.dir/field.f.o

CMakeFiles/CADD.dir/mod_crack.f.o: CMakeFiles/CADD.dir/flags.make
CMakeFiles/CADD.dir/mod_crack.f.o: ../mod_crack.f
	$(CMAKE_COMMAND) -E cmake_progress_report /home/srinath/CADD_all/CADD_ben2/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building Fortran object CMakeFiles/CADD.dir/mod_crack.f.o"
	/opt/intel/bin/ifort  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /home/srinath/CADD_all/CADD_ben2/mod_crack.f -o CMakeFiles/CADD.dir/mod_crack.f.o

CMakeFiles/CADD.dir/mod_crack.f.o.requires:
.PHONY : CMakeFiles/CADD.dir/mod_crack.f.o.requires

CMakeFiles/CADD.dir/mod_crack.f.o.provides: CMakeFiles/CADD.dir/mod_crack.f.o.requires
	$(MAKE) -f CMakeFiles/CADD.dir/build.make CMakeFiles/CADD.dir/mod_crack.f.o.provides.build
.PHONY : CMakeFiles/CADD.dir/mod_crack.f.o.provides

CMakeFiles/CADD.dir/mod_crack.f.o.provides.build: CMakeFiles/CADD.dir/mod_crack.f.o

CMakeFiles/CADD.dir/mesh.f.o: CMakeFiles/CADD.dir/flags.make
CMakeFiles/CADD.dir/mesh.f.o: ../mesh.f
	$(CMAKE_COMMAND) -E cmake_progress_report /home/srinath/CADD_all/CADD_ben2/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building Fortran object CMakeFiles/CADD.dir/mesh.f.o"
	/opt/intel/bin/ifort  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /home/srinath/CADD_all/CADD_ben2/mesh.f -o CMakeFiles/CADD.dir/mesh.f.o

CMakeFiles/CADD.dir/mesh.f.o.requires:
.PHONY : CMakeFiles/CADD.dir/mesh.f.o.requires

CMakeFiles/CADD.dir/mesh.f.o.provides: CMakeFiles/CADD.dir/mesh.f.o.requires
	$(MAKE) -f CMakeFiles/CADD.dir/build.make CMakeFiles/CADD.dir/mesh.f.o.provides.build
.PHONY : CMakeFiles/CADD.dir/mesh.f.o.provides

CMakeFiles/CADD.dir/mesh.f.o.provides.build: CMakeFiles/CADD.dir/mesh.f.o

CMakeFiles/CADD.dir/move_crack.f.o: CMakeFiles/CADD.dir/flags.make
CMakeFiles/CADD.dir/move_crack.f.o: ../move_crack.f
	$(CMAKE_COMMAND) -E cmake_progress_report /home/srinath/CADD_all/CADD_ben2/build/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building Fortran object CMakeFiles/CADD.dir/move_crack.f.o"
	/opt/intel/bin/ifort  $(Fortran_DEFINES) $(Fortran_FLAGS) -c /home/srinath/CADD_all/CADD_ben2/move_crack.f -o CMakeFiles/CADD.dir/move_crack.f.o

CMakeFiles/CADD.dir/move_crack.f.o.requires:
.PHONY : CMakeFiles/CADD.dir/move_crack.f.o.requires

CMakeFiles/CADD.dir/move_crack.f.o.provides: CMakeFiles/CADD.dir/move_crack.f.o.requires
	$(MAKE) -f CMakeFiles/CADD.dir/build.make CMakeFiles/CADD.dir/move_crack.f.o.provides.build
.PHONY : CMakeFiles/CADD.dir/move_crack.f.o.provides

CMakeFiles/CADD.dir/move_crack.f.o.provides.build: CMakeFiles/CADD.dir/move_crack.f.o

# Object files for target CADD
CADD_OBJECTS = \
"CMakeFiles/CADD.dir/field.f.o" \
"CMakeFiles/CADD.dir/mod_crack.f.o" \
"CMakeFiles/CADD.dir/mesh.f.o" \
"CMakeFiles/CADD.dir/move_crack.f.o"

# External object files for target CADD
CADD_EXTERNAL_OBJECTS =

CADD: CMakeFiles/CADD.dir/field.f.o
CADD: CMakeFiles/CADD.dir/mod_crack.f.o
CADD: CMakeFiles/CADD.dir/mesh.f.o
CADD: CMakeFiles/CADD.dir/move_crack.f.o
CADD: CMakeFiles/CADD.dir/build.make
CADD: Disl/libdislo.so
CADD: Modular/libmdlro.so
CADD: CMakeFiles/CADD.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking Fortran executable CADD"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/CADD.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/CADD.dir/build: CADD
.PHONY : CMakeFiles/CADD.dir/build

CMakeFiles/CADD.dir/requires: CMakeFiles/CADD.dir/field.f.o.requires
CMakeFiles/CADD.dir/requires: CMakeFiles/CADD.dir/mod_crack.f.o.requires
CMakeFiles/CADD.dir/requires: CMakeFiles/CADD.dir/mesh.f.o.requires
CMakeFiles/CADD.dir/requires: CMakeFiles/CADD.dir/move_crack.f.o.requires
.PHONY : CMakeFiles/CADD.dir/requires

CMakeFiles/CADD.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/CADD.dir/cmake_clean.cmake
.PHONY : CMakeFiles/CADD.dir/clean

CMakeFiles/CADD.dir/depend:
	cd /home/srinath/CADD_all/CADD_ben2/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/srinath/CADD_all/CADD_ben2 /home/srinath/CADD_all/CADD_ben2 /home/srinath/CADD_all/CADD_ben2/build /home/srinath/CADD_all/CADD_ben2/build /home/srinath/CADD_all/CADD_ben2/build/CMakeFiles/CADD.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/CADD.dir/depend
