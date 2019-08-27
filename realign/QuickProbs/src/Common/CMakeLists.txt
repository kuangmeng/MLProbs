cmake_minimum_required (VERSION 2.6)
project (Common)

set(SRC_FILES 
	Array.h
	dbgnew.h
	deterministic_random.cpp
	deterministic_random.h
	dirent.h
	mathex.cpp
	mathex.h
	Log.h
	Log.cpp
	Named.h
	PlacementAllocator.h
	Printable.h
	ProgramOptions.cpp
	ProgramOptions.h
	StatisticsProvider.cpp
	StatisticsProvider.h
	Timer.h
	MemoryTools.h
	WindowsMemoryTools.h
	WindowsMemoryTools.cpp
	UnixMemoryTools.h
	UnixMemoryTools.cpp)

source_group("" FILES ${SRC_FILES})
add_library(Common STATIC ${SRC_FILES})