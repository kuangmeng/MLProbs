cmake_minimum_required (VERSION 2.6)
project (tcmalloc)

set(HEADER_FILES
	addressmap-inl.h
    base/atomicops-internals-linuxppc.h
    base/atomicops-internals-macosx.h
    base/atomicops-internals-x86.h
    base/atomicops.h
    base/basictypes.h
    base/commandlineflags.h
    base/googleinit.h
    base/logging.h
    base/low_level_alloc.h
    base/spinlock.h
    base/spinlock_internal.h
    base/spinlock_linux-inl.h
    base/spinlock_posix-inl.h
    base/spinlock_win32-inl.h
    base/stl_allocator.h
    base/sysinfo.h
    base/thread_annotations.h
    central_freelist.h
    common.h
    gperftools/heap-checker.h
    gperftools/heap-profiler.h
    gperftools/malloc_extension.h
    gperftools/malloc_hook.h
    gperftools/profiler.h
    gperftools/stacktrace.h
    heap-profile-table.h
    internal_logging.h
    malloc_hook-inl.h
    memory_region_map.h
    packed-cache-inl.h
    pagemap.h
    page_heap.h
    page_heap_allocator.h
    raw_printer.h
    sampler.h
    span.h
    stacktrace_win32-inl.h
    stack_trace_table.h
    static_vars.h
    symbolize.h
    system-alloc.h
    tcmalloc.h
    thread_cache.h
    windows/config.h
    windows/mini_disassembler.h
    windows/mini_disassembler_types.h
    windows/port.h
    windows/preamble_patcher.h
)

set(SOURCE_FILES 
	central_freelist.cc
	base/dynamic_annotations.c
	heap-profile-table.cc
	symbolize.cc
	windows/ia32_modrm_map.cc
	windows/ia32_opcode_map.cc
	common.cc
	internal_logging.cc
	base/logging.cc
	base/low_level_alloc.cc
	malloc_extension.cc
	malloc_hook.cc
	memory_region_map.cc
	windows/mini_disassembler.cc
	page_heap.cc
	sampler.cc
	windows/patch_functions.cc
	windows/port.cc
	windows/preamble_patcher.cc
	windows/preamble_patcher_with_stub.cc
	windows/system-alloc.cc
	raw_printer.cc
	span.cc
	stacktrace.cc
	stack_trace_table.cc
	static_vars.cc
	base/spinlock.cc
	base/spinlock_internal.cc
	base/sysinfo.cc
	thread_cache.cc
	fake_stacktrace_scope.cc
)

include_directories(".")
include_directories("./windows")

source_group("Header Files" FILES ${HEADER_FILES})
source_group("Source Files" FILES ${SOURCE_FILES})

add_library(tcmalloc ${HEADER_FILES} ${SOURCE_FILES})


#if (WIN32)
#	add_custom_command(TARGET KernelRepository 
#	    PRE_BUILD
#	COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR}
#	    COMMAND python process_kernel_files.py
#	)
#endif()
