cmake_minimum_required (VERSION 2.6)
project (KernelAlignment)

set(MAIN_FILES
	KernelFactory.cpp
	KernelFactory.h
	KernelRepository.cpp
	KernelRepository.h
	process_kernel_files.py
)

set(KERNEL_FILES 
	Kernels/AuxiliaryStructures.cl
	Kernels/BinarySearch.cl
	Kernels/Finalization.cl
	Kernels/Finalization_long.cl
	Kernels/LogarithmUtils.cl
	Kernels/MemoryTransfers.cl
	Kernels/Partition.cl
	Kernels/PartitionLogarithm.cl
	Kernels/Partition_long.cl
	Kernels/PartitionLogarithm_long.cl
	Kernels/Probabilistic.cl
	Kernels/Probabilistic_long.cl
	Kernels/Random.cl
	Kernels/ScoreType.cl
	Kernels/SparseMatrix.cl
	Kernels/SparseMatrix_long.cl
	Kernels/Utils.cl

	Kernels/MultipleConsistency_old.cl
	Kernels/MultipleConsistencyUtils.cl
	Kernels/MultiplePartition.cl
	Kernels/MultipleProbabilistic.cl
	Kernels/MultipleFinalization.cl
	Kernels/MultipleSparse.cl
)

set(MSA_FILES
	Multiple/ITasksWave.h
	Multiple/KernelMSA.cpp
	Multiple/KernelMSA.h
	Multiple/KernelSet.h
	Multiple/KernelSet.cpp
	Multiple/PosteriorTasksWave.cpp
	Multiple/PosteriorTasksWave.h
	Multiple/QuickPosteriorStage.cpp
	Multiple/QuickPosteriorStage.h
	Multiple/QuickConsistencyStage.cpp
	Multiple/QuickConsistencyStage.h
	Multiple/RelaxationKernelWrapper.h
	Multiple/RelaxationSector.cpp
	Multiple/RelaxationSector.h
	
	)
	
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${Extra_INCLUDE_DIRS})

source_group("" FILES ${MAIN_FILES})
source_group("Multiple" FILES ${MSA_FILES})
source_group("Kernels" FILES ${KERNEL_FILES})

add_library(KernelAlignment STATIC ${MAIN_FILES} ${KERNEL_FILES} ${MSA_FILES})
target_link_libraries(KernelAlignment 
	debug ${TCMalloc_LIBRARIES_DEBUG} 
	optimized ${TCMalloc_LIBRARIES_RELEASE}
	general ${OpenCL_LIBRARIES} Alignment Common Hardware)

#if (WIN32)
#	add_custom_command(TARGET KernelRepository 
#	    PRE_BUILD
#	COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR}
#	    COMMAND python process_kernel_files.py
#	)
#endif()
