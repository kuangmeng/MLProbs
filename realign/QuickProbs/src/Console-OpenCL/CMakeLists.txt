cmake_minimum_required (VERSION 2.6)
project (Console-OpenCl)

set(SRC_FILES main-opencl.cpp)

set (EXE_NAME "quickprobs-cl")
if (WIN32)
	set (EXE_NAME "${EXE_NAME}-win32")
elseif(OSX)
	set (EXE_NAME "${EXE_NAME}-osx")
endif()

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${Extra_INCLUDE_DIRS})

source_group("" FILES ${SRC_FILES})
add_executable(Console-OpenCl ${SRC_FILES})
target_link_libraries(Console-OpenCl 
	debug ${TCMalloc_LIBRARIES_DEBUG} 
	optimized ${TCMalloc_LIBRARIES_RELEASE}
	general ${System_LIBRARIES} ${OpenCL_LIBRARIES} Alignment KernelAlignment Common Hardware)

set_target_properties(Console-OpenCl PROPERTIES OUTPUT_NAME ${EXE_NAME} ) 
set(EXECUTABLE_OUTPUT_PATH ${OUT_PATH})

