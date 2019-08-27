cmake_minimum_required (VERSION 2.6)
project (Hardware)

set(SRC_FILES 
	EventHandler.h
	DeviceInfo.cpp
	DeviceInfo.h
	DeviceWrapper.h
	DeviceWrapper.cpp
	OpenCl.cpp
	OpenCl.h
	OpenCLError.h
	Buffer.cpp
	Buffer.h
	Kernel.h 
	PlatformInfo.h)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${Extra_INCLUDE_DIRS})	
			
source_group("" FILES ${SRC_FILES})
add_library(Hardware STATIC ${SRC_FILES})
target_link_libraries(Hardware ${OpenCL_LIBRARIES} Common)
