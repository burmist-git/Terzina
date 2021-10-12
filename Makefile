# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

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
CMAKE_SOURCE_DIR = /home/burmist/geant4_workdir/lhcBrich

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/burmist/geant4_workdir/lhcBrich

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running interactive CMake command-line interface..."
	/usr/bin/cmake -i .
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target install
install: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/usr/bin/cmake -P cmake_install.cmake
.PHONY : install

# Special rule for the target install
install/fast: preinstall/fast
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Install the project..."
	/usr/bin/cmake -P cmake_install.cmake
.PHONY : install/fast

# Special rule for the target install/local
install/local: preinstall
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Installing only the local directory..."
	/usr/bin/cmake -DCMAKE_INSTALL_LOCAL_ONLY=1 -P cmake_install.cmake
.PHONY : install/local

# Special rule for the target install/local
install/local/fast: install/local
.PHONY : install/local/fast

# Special rule for the target list_install_components
list_install_components:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Available install components are: \"Unspecified\""
.PHONY : list_install_components

# Special rule for the target list_install_components
list_install_components/fast: list_install_components
.PHONY : list_install_components/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/burmist/geant4_workdir/lhcBrich/CMakeFiles /home/burmist/geant4_workdir/lhcBrich/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/burmist/geant4_workdir/lhcBrich/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named lhcBrich

# Build rule for target.
lhcBrich: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 lhcBrich
.PHONY : lhcBrich

# fast build rule for target.
lhcBrich/fast:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/build
.PHONY : lhcBrich/fast

# Manual pre-install relink rule for target.
lhcBrich/preinstall:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/preinstall
.PHONY : lhcBrich/preinstall

lhcBrich.o: lhcBrich.cc.o
.PHONY : lhcBrich.o

# target to build an object file
lhcBrich.cc.o:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/lhcBrich.cc.o
.PHONY : lhcBrich.cc.o

lhcBrich.i: lhcBrich.cc.i
.PHONY : lhcBrich.i

# target to preprocess a source file
lhcBrich.cc.i:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/lhcBrich.cc.i
.PHONY : lhcBrich.cc.i

lhcBrich.s: lhcBrich.cc.s
.PHONY : lhcBrich.s

# target to generate assembly for a file
lhcBrich.cc.s:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/lhcBrich.cc.s
.PHONY : lhcBrich.cc.s

src/DetectorConstruction.o: src/DetectorConstruction.cc.o
.PHONY : src/DetectorConstruction.o

# target to build an object file
src/DetectorConstruction.cc.o:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/DetectorConstruction.cc.o
.PHONY : src/DetectorConstruction.cc.o

src/DetectorConstruction.i: src/DetectorConstruction.cc.i
.PHONY : src/DetectorConstruction.i

# target to preprocess a source file
src/DetectorConstruction.cc.i:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/DetectorConstruction.cc.i
.PHONY : src/DetectorConstruction.cc.i

src/DetectorConstruction.s: src/DetectorConstruction.cc.s
.PHONY : src/DetectorConstruction.s

# target to generate assembly for a file
src/DetectorConstruction.cc.s:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/DetectorConstruction.cc.s
.PHONY : src/DetectorConstruction.cc.s

src/EventAction.o: src/EventAction.cc.o
.PHONY : src/EventAction.o

# target to build an object file
src/EventAction.cc.o:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/EventAction.cc.o
.PHONY : src/EventAction.cc.o

src/EventAction.i: src/EventAction.cc.i
.PHONY : src/EventAction.i

# target to preprocess a source file
src/EventAction.cc.i:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/EventAction.cc.i
.PHONY : src/EventAction.cc.i

src/EventAction.s: src/EventAction.cc.s
.PHONY : src/EventAction.s

# target to generate assembly for a file
src/EventAction.cc.s:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/EventAction.cc.s
.PHONY : src/EventAction.cc.s

src/HitMy.o: src/HitMy.cc.o
.PHONY : src/HitMy.o

# target to build an object file
src/HitMy.cc.o:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/HitMy.cc.o
.PHONY : src/HitMy.cc.o

src/HitMy.i: src/HitMy.cc.i
.PHONY : src/HitMy.i

# target to preprocess a source file
src/HitMy.cc.i:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/HitMy.cc.i
.PHONY : src/HitMy.cc.i

src/HitMy.s: src/HitMy.cc.s
.PHONY : src/HitMy.s

# target to generate assembly for a file
src/HitMy.cc.s:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/HitMy.cc.s
.PHONY : src/HitMy.cc.s

src/PhysicsList.o: src/PhysicsList.cc.o
.PHONY : src/PhysicsList.o

# target to build an object file
src/PhysicsList.cc.o:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/PhysicsList.cc.o
.PHONY : src/PhysicsList.cc.o

src/PhysicsList.i: src/PhysicsList.cc.i
.PHONY : src/PhysicsList.i

# target to preprocess a source file
src/PhysicsList.cc.i:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/PhysicsList.cc.i
.PHONY : src/PhysicsList.cc.i

src/PhysicsList.s: src/PhysicsList.cc.s
.PHONY : src/PhysicsList.s

# target to generate assembly for a file
src/PhysicsList.cc.s:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/PhysicsList.cc.s
.PHONY : src/PhysicsList.cc.s

src/PrimaryGeneratorAction.o: src/PrimaryGeneratorAction.cc.o
.PHONY : src/PrimaryGeneratorAction.o

# target to build an object file
src/PrimaryGeneratorAction.cc.o:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/PrimaryGeneratorAction.cc.o
.PHONY : src/PrimaryGeneratorAction.cc.o

src/PrimaryGeneratorAction.i: src/PrimaryGeneratorAction.cc.i
.PHONY : src/PrimaryGeneratorAction.i

# target to preprocess a source file
src/PrimaryGeneratorAction.cc.i:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/PrimaryGeneratorAction.cc.i
.PHONY : src/PrimaryGeneratorAction.cc.i

src/PrimaryGeneratorAction.s: src/PrimaryGeneratorAction.cc.s
.PHONY : src/PrimaryGeneratorAction.s

# target to generate assembly for a file
src/PrimaryGeneratorAction.cc.s:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/PrimaryGeneratorAction.cc.s
.PHONY : src/PrimaryGeneratorAction.cc.s

src/RunAction.o: src/RunAction.cc.o
.PHONY : src/RunAction.o

# target to build an object file
src/RunAction.cc.o:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/RunAction.cc.o
.PHONY : src/RunAction.cc.o

src/RunAction.i: src/RunAction.cc.i
.PHONY : src/RunAction.i

# target to preprocess a source file
src/RunAction.cc.i:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/RunAction.cc.i
.PHONY : src/RunAction.cc.i

src/RunAction.s: src/RunAction.cc.s
.PHONY : src/RunAction.s

# target to generate assembly for a file
src/RunAction.cc.s:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/RunAction.cc.s
.PHONY : src/RunAction.cc.s

src/SensitiveDetector.o: src/SensitiveDetector.cc.o
.PHONY : src/SensitiveDetector.o

# target to build an object file
src/SensitiveDetector.cc.o:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/SensitiveDetector.cc.o
.PHONY : src/SensitiveDetector.cc.o

src/SensitiveDetector.i: src/SensitiveDetector.cc.i
.PHONY : src/SensitiveDetector.i

# target to preprocess a source file
src/SensitiveDetector.cc.i:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/SensitiveDetector.cc.i
.PHONY : src/SensitiveDetector.cc.i

src/SensitiveDetector.s: src/SensitiveDetector.cc.s
.PHONY : src/SensitiveDetector.s

# target to generate assembly for a file
src/SensitiveDetector.cc.s:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/SensitiveDetector.cc.s
.PHONY : src/SensitiveDetector.cc.s

src/StackingAction.o: src/StackingAction.cc.o
.PHONY : src/StackingAction.o

# target to build an object file
src/StackingAction.cc.o:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/StackingAction.cc.o
.PHONY : src/StackingAction.cc.o

src/StackingAction.i: src/StackingAction.cc.i
.PHONY : src/StackingAction.i

# target to preprocess a source file
src/StackingAction.cc.i:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/StackingAction.cc.i
.PHONY : src/StackingAction.cc.i

src/StackingAction.s: src/StackingAction.cc.s
.PHONY : src/StackingAction.s

# target to generate assembly for a file
src/StackingAction.cc.s:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/StackingAction.cc.s
.PHONY : src/StackingAction.cc.s

src/SteppingAction.o: src/SteppingAction.cc.o
.PHONY : src/SteppingAction.o

# target to build an object file
src/SteppingAction.cc.o:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/SteppingAction.cc.o
.PHONY : src/SteppingAction.cc.o

src/SteppingAction.i: src/SteppingAction.cc.i
.PHONY : src/SteppingAction.i

# target to preprocess a source file
src/SteppingAction.cc.i:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/SteppingAction.cc.i
.PHONY : src/SteppingAction.cc.i

src/SteppingAction.s: src/SteppingAction.cc.s
.PHONY : src/SteppingAction.s

# target to generate assembly for a file
src/SteppingAction.cc.s:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/SteppingAction.cc.s
.PHONY : src/SteppingAction.cc.s

src/SteppingVerbose.o: src/SteppingVerbose.cc.o
.PHONY : src/SteppingVerbose.o

# target to build an object file
src/SteppingVerbose.cc.o:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/SteppingVerbose.cc.o
.PHONY : src/SteppingVerbose.cc.o

src/SteppingVerbose.i: src/SteppingVerbose.cc.i
.PHONY : src/SteppingVerbose.i

# target to preprocess a source file
src/SteppingVerbose.cc.i:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/SteppingVerbose.cc.i
.PHONY : src/SteppingVerbose.cc.i

src/SteppingVerbose.s: src/SteppingVerbose.cc.s
.PHONY : src/SteppingVerbose.s

# target to generate assembly for a file
src/SteppingVerbose.cc.s:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/SteppingVerbose.cc.s
.PHONY : src/SteppingVerbose.cc.s

src/TrackInformation.o: src/TrackInformation.cc.o
.PHONY : src/TrackInformation.o

# target to build an object file
src/TrackInformation.cc.o:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/TrackInformation.cc.o
.PHONY : src/TrackInformation.cc.o

src/TrackInformation.i: src/TrackInformation.cc.i
.PHONY : src/TrackInformation.i

# target to preprocess a source file
src/TrackInformation.cc.i:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/TrackInformation.cc.i
.PHONY : src/TrackInformation.cc.i

src/TrackInformation.s: src/TrackInformation.cc.s
.PHONY : src/TrackInformation.s

# target to generate assembly for a file
src/TrackInformation.cc.s:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/TrackInformation.cc.s
.PHONY : src/TrackInformation.cc.s

src/TrackingAction.o: src/TrackingAction.cc.o
.PHONY : src/TrackingAction.o

# target to build an object file
src/TrackingAction.cc.o:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/TrackingAction.cc.o
.PHONY : src/TrackingAction.cc.o

src/TrackingAction.i: src/TrackingAction.cc.i
.PHONY : src/TrackingAction.i

# target to preprocess a source file
src/TrackingAction.cc.i:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/TrackingAction.cc.i
.PHONY : src/TrackingAction.cc.i

src/TrackingAction.s: src/TrackingAction.cc.s
.PHONY : src/TrackingAction.s

# target to generate assembly for a file
src/TrackingAction.cc.s:
	$(MAKE) -f CMakeFiles/lhcBrich.dir/build.make CMakeFiles/lhcBrich.dir/src/TrackingAction.cc.s
.PHONY : src/TrackingAction.cc.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... edit_cache"
	@echo "... install"
	@echo "... install/local"
	@echo "... lhcBrich"
	@echo "... list_install_components"
	@echo "... rebuild_cache"
	@echo "... lhcBrich.o"
	@echo "... lhcBrich.i"
	@echo "... lhcBrich.s"
	@echo "... src/DetectorConstruction.o"
	@echo "... src/DetectorConstruction.i"
	@echo "... src/DetectorConstruction.s"
	@echo "... src/EventAction.o"
	@echo "... src/EventAction.i"
	@echo "... src/EventAction.s"
	@echo "... src/HitMy.o"
	@echo "... src/HitMy.i"
	@echo "... src/HitMy.s"
	@echo "... src/PhysicsList.o"
	@echo "... src/PhysicsList.i"
	@echo "... src/PhysicsList.s"
	@echo "... src/PrimaryGeneratorAction.o"
	@echo "... src/PrimaryGeneratorAction.i"
	@echo "... src/PrimaryGeneratorAction.s"
	@echo "... src/RunAction.o"
	@echo "... src/RunAction.i"
	@echo "... src/RunAction.s"
	@echo "... src/SensitiveDetector.o"
	@echo "... src/SensitiveDetector.i"
	@echo "... src/SensitiveDetector.s"
	@echo "... src/StackingAction.o"
	@echo "... src/StackingAction.i"
	@echo "... src/StackingAction.s"
	@echo "... src/SteppingAction.o"
	@echo "... src/SteppingAction.i"
	@echo "... src/SteppingAction.s"
	@echo "... src/SteppingVerbose.o"
	@echo "... src/SteppingVerbose.i"
	@echo "... src/SteppingVerbose.s"
	@echo "... src/TrackInformation.o"
	@echo "... src/TrackInformation.i"
	@echo "... src/TrackInformation.s"
	@echo "... src/TrackingAction.o"
	@echo "... src/TrackingAction.i"
	@echo "... src/TrackingAction.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system
