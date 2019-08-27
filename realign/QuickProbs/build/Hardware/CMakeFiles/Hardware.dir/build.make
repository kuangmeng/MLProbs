# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_SOURCE_DIR = /home/lg28/桌面/QuickProbs/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/lg28/桌面/QuickProbs/build

# Include any dependencies generated for this target.
include Hardware/CMakeFiles/Hardware.dir/depend.make

# Include the progress variables for this target.
include Hardware/CMakeFiles/Hardware.dir/progress.make

# Include the compile flags for this target's objects.
include Hardware/CMakeFiles/Hardware.dir/flags.make

Hardware/CMakeFiles/Hardware.dir/DeviceInfo.cpp.o: Hardware/CMakeFiles/Hardware.dir/flags.make
Hardware/CMakeFiles/Hardware.dir/DeviceInfo.cpp.o: /home/lg28/桌面/QuickProbs/src/Hardware/DeviceInfo.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lg28/桌面/QuickProbs/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object Hardware/CMakeFiles/Hardware.dir/DeviceInfo.cpp.o"
	cd /home/lg28/桌面/QuickProbs/build/Hardware && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Hardware.dir/DeviceInfo.cpp.o -c /home/lg28/桌面/QuickProbs/src/Hardware/DeviceInfo.cpp

Hardware/CMakeFiles/Hardware.dir/DeviceInfo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Hardware.dir/DeviceInfo.cpp.i"
	cd /home/lg28/桌面/QuickProbs/build/Hardware && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lg28/桌面/QuickProbs/src/Hardware/DeviceInfo.cpp > CMakeFiles/Hardware.dir/DeviceInfo.cpp.i

Hardware/CMakeFiles/Hardware.dir/DeviceInfo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Hardware.dir/DeviceInfo.cpp.s"
	cd /home/lg28/桌面/QuickProbs/build/Hardware && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lg28/桌面/QuickProbs/src/Hardware/DeviceInfo.cpp -o CMakeFiles/Hardware.dir/DeviceInfo.cpp.s

Hardware/CMakeFiles/Hardware.dir/DeviceInfo.cpp.o.requires:

.PHONY : Hardware/CMakeFiles/Hardware.dir/DeviceInfo.cpp.o.requires

Hardware/CMakeFiles/Hardware.dir/DeviceInfo.cpp.o.provides: Hardware/CMakeFiles/Hardware.dir/DeviceInfo.cpp.o.requires
	$(MAKE) -f Hardware/CMakeFiles/Hardware.dir/build.make Hardware/CMakeFiles/Hardware.dir/DeviceInfo.cpp.o.provides.build
.PHONY : Hardware/CMakeFiles/Hardware.dir/DeviceInfo.cpp.o.provides

Hardware/CMakeFiles/Hardware.dir/DeviceInfo.cpp.o.provides.build: Hardware/CMakeFiles/Hardware.dir/DeviceInfo.cpp.o


Hardware/CMakeFiles/Hardware.dir/DeviceWrapper.cpp.o: Hardware/CMakeFiles/Hardware.dir/flags.make
Hardware/CMakeFiles/Hardware.dir/DeviceWrapper.cpp.o: /home/lg28/桌面/QuickProbs/src/Hardware/DeviceWrapper.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lg28/桌面/QuickProbs/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object Hardware/CMakeFiles/Hardware.dir/DeviceWrapper.cpp.o"
	cd /home/lg28/桌面/QuickProbs/build/Hardware && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Hardware.dir/DeviceWrapper.cpp.o -c /home/lg28/桌面/QuickProbs/src/Hardware/DeviceWrapper.cpp

Hardware/CMakeFiles/Hardware.dir/DeviceWrapper.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Hardware.dir/DeviceWrapper.cpp.i"
	cd /home/lg28/桌面/QuickProbs/build/Hardware && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lg28/桌面/QuickProbs/src/Hardware/DeviceWrapper.cpp > CMakeFiles/Hardware.dir/DeviceWrapper.cpp.i

Hardware/CMakeFiles/Hardware.dir/DeviceWrapper.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Hardware.dir/DeviceWrapper.cpp.s"
	cd /home/lg28/桌面/QuickProbs/build/Hardware && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lg28/桌面/QuickProbs/src/Hardware/DeviceWrapper.cpp -o CMakeFiles/Hardware.dir/DeviceWrapper.cpp.s

Hardware/CMakeFiles/Hardware.dir/DeviceWrapper.cpp.o.requires:

.PHONY : Hardware/CMakeFiles/Hardware.dir/DeviceWrapper.cpp.o.requires

Hardware/CMakeFiles/Hardware.dir/DeviceWrapper.cpp.o.provides: Hardware/CMakeFiles/Hardware.dir/DeviceWrapper.cpp.o.requires
	$(MAKE) -f Hardware/CMakeFiles/Hardware.dir/build.make Hardware/CMakeFiles/Hardware.dir/DeviceWrapper.cpp.o.provides.build
.PHONY : Hardware/CMakeFiles/Hardware.dir/DeviceWrapper.cpp.o.provides

Hardware/CMakeFiles/Hardware.dir/DeviceWrapper.cpp.o.provides.build: Hardware/CMakeFiles/Hardware.dir/DeviceWrapper.cpp.o


Hardware/CMakeFiles/Hardware.dir/OpenCl.cpp.o: Hardware/CMakeFiles/Hardware.dir/flags.make
Hardware/CMakeFiles/Hardware.dir/OpenCl.cpp.o: /home/lg28/桌面/QuickProbs/src/Hardware/OpenCl.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lg28/桌面/QuickProbs/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object Hardware/CMakeFiles/Hardware.dir/OpenCl.cpp.o"
	cd /home/lg28/桌面/QuickProbs/build/Hardware && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Hardware.dir/OpenCl.cpp.o -c /home/lg28/桌面/QuickProbs/src/Hardware/OpenCl.cpp

Hardware/CMakeFiles/Hardware.dir/OpenCl.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Hardware.dir/OpenCl.cpp.i"
	cd /home/lg28/桌面/QuickProbs/build/Hardware && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lg28/桌面/QuickProbs/src/Hardware/OpenCl.cpp > CMakeFiles/Hardware.dir/OpenCl.cpp.i

Hardware/CMakeFiles/Hardware.dir/OpenCl.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Hardware.dir/OpenCl.cpp.s"
	cd /home/lg28/桌面/QuickProbs/build/Hardware && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lg28/桌面/QuickProbs/src/Hardware/OpenCl.cpp -o CMakeFiles/Hardware.dir/OpenCl.cpp.s

Hardware/CMakeFiles/Hardware.dir/OpenCl.cpp.o.requires:

.PHONY : Hardware/CMakeFiles/Hardware.dir/OpenCl.cpp.o.requires

Hardware/CMakeFiles/Hardware.dir/OpenCl.cpp.o.provides: Hardware/CMakeFiles/Hardware.dir/OpenCl.cpp.o.requires
	$(MAKE) -f Hardware/CMakeFiles/Hardware.dir/build.make Hardware/CMakeFiles/Hardware.dir/OpenCl.cpp.o.provides.build
.PHONY : Hardware/CMakeFiles/Hardware.dir/OpenCl.cpp.o.provides

Hardware/CMakeFiles/Hardware.dir/OpenCl.cpp.o.provides.build: Hardware/CMakeFiles/Hardware.dir/OpenCl.cpp.o


Hardware/CMakeFiles/Hardware.dir/Buffer.cpp.o: Hardware/CMakeFiles/Hardware.dir/flags.make
Hardware/CMakeFiles/Hardware.dir/Buffer.cpp.o: /home/lg28/桌面/QuickProbs/src/Hardware/Buffer.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/lg28/桌面/QuickProbs/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object Hardware/CMakeFiles/Hardware.dir/Buffer.cpp.o"
	cd /home/lg28/桌面/QuickProbs/build/Hardware && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/Hardware.dir/Buffer.cpp.o -c /home/lg28/桌面/QuickProbs/src/Hardware/Buffer.cpp

Hardware/CMakeFiles/Hardware.dir/Buffer.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Hardware.dir/Buffer.cpp.i"
	cd /home/lg28/桌面/QuickProbs/build/Hardware && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/lg28/桌面/QuickProbs/src/Hardware/Buffer.cpp > CMakeFiles/Hardware.dir/Buffer.cpp.i

Hardware/CMakeFiles/Hardware.dir/Buffer.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Hardware.dir/Buffer.cpp.s"
	cd /home/lg28/桌面/QuickProbs/build/Hardware && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/lg28/桌面/QuickProbs/src/Hardware/Buffer.cpp -o CMakeFiles/Hardware.dir/Buffer.cpp.s

Hardware/CMakeFiles/Hardware.dir/Buffer.cpp.o.requires:

.PHONY : Hardware/CMakeFiles/Hardware.dir/Buffer.cpp.o.requires

Hardware/CMakeFiles/Hardware.dir/Buffer.cpp.o.provides: Hardware/CMakeFiles/Hardware.dir/Buffer.cpp.o.requires
	$(MAKE) -f Hardware/CMakeFiles/Hardware.dir/build.make Hardware/CMakeFiles/Hardware.dir/Buffer.cpp.o.provides.build
.PHONY : Hardware/CMakeFiles/Hardware.dir/Buffer.cpp.o.provides

Hardware/CMakeFiles/Hardware.dir/Buffer.cpp.o.provides.build: Hardware/CMakeFiles/Hardware.dir/Buffer.cpp.o


# Object files for target Hardware
Hardware_OBJECTS = \
"CMakeFiles/Hardware.dir/DeviceInfo.cpp.o" \
"CMakeFiles/Hardware.dir/DeviceWrapper.cpp.o" \
"CMakeFiles/Hardware.dir/OpenCl.cpp.o" \
"CMakeFiles/Hardware.dir/Buffer.cpp.o"

# External object files for target Hardware
Hardware_EXTERNAL_OBJECTS =

Hardware/libHardware.a: Hardware/CMakeFiles/Hardware.dir/DeviceInfo.cpp.o
Hardware/libHardware.a: Hardware/CMakeFiles/Hardware.dir/DeviceWrapper.cpp.o
Hardware/libHardware.a: Hardware/CMakeFiles/Hardware.dir/OpenCl.cpp.o
Hardware/libHardware.a: Hardware/CMakeFiles/Hardware.dir/Buffer.cpp.o
Hardware/libHardware.a: Hardware/CMakeFiles/Hardware.dir/build.make
Hardware/libHardware.a: Hardware/CMakeFiles/Hardware.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/lg28/桌面/QuickProbs/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX static library libHardware.a"
	cd /home/lg28/桌面/QuickProbs/build/Hardware && $(CMAKE_COMMAND) -P CMakeFiles/Hardware.dir/cmake_clean_target.cmake
	cd /home/lg28/桌面/QuickProbs/build/Hardware && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/Hardware.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Hardware/CMakeFiles/Hardware.dir/build: Hardware/libHardware.a

.PHONY : Hardware/CMakeFiles/Hardware.dir/build

Hardware/CMakeFiles/Hardware.dir/requires: Hardware/CMakeFiles/Hardware.dir/DeviceInfo.cpp.o.requires
Hardware/CMakeFiles/Hardware.dir/requires: Hardware/CMakeFiles/Hardware.dir/DeviceWrapper.cpp.o.requires
Hardware/CMakeFiles/Hardware.dir/requires: Hardware/CMakeFiles/Hardware.dir/OpenCl.cpp.o.requires
Hardware/CMakeFiles/Hardware.dir/requires: Hardware/CMakeFiles/Hardware.dir/Buffer.cpp.o.requires

.PHONY : Hardware/CMakeFiles/Hardware.dir/requires

Hardware/CMakeFiles/Hardware.dir/clean:
	cd /home/lg28/桌面/QuickProbs/build/Hardware && $(CMAKE_COMMAND) -P CMakeFiles/Hardware.dir/cmake_clean.cmake
.PHONY : Hardware/CMakeFiles/Hardware.dir/clean

Hardware/CMakeFiles/Hardware.dir/depend:
	cd /home/lg28/桌面/QuickProbs/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/lg28/桌面/QuickProbs/src /home/lg28/桌面/QuickProbs/src/Hardware /home/lg28/桌面/QuickProbs/build /home/lg28/桌面/QuickProbs/build/Hardware /home/lg28/桌面/QuickProbs/build/Hardware/CMakeFiles/Hardware.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Hardware/CMakeFiles/Hardware.dir/depend
