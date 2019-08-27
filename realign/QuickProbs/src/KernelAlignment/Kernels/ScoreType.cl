#if defined(cl_khr_fp64) //&& __OPENCL_VERSION__ == 10	// Khronos extension available?
	#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64) //&& __OPENCL_VERSION__ == 10 // AMD extension available?
	#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

#define assert(x)	 // empty function

// some constants
#define NumMatchStates		1
#define NumInsertStates 	2
#define NumMatrixTypes		5
#define SymbolCountSmall	26
#define depth				NumMatrixTypes

// define special memory qualifiers
#ifdef COPY_PROBABILISTIC_PARAMS
	#define HMM_MEM local
	#define MAX_SIZE(x)  //__attribute__((max_constant_size (x)))
#else
	#define HMM_MEM global
	#define MAX_SIZE(x) 
#endif

#ifdef COPY_PARTITION_PARAMS
	#define PART_MEM local
	#define MAX_SIZE(x)  //__attribute__((max_constant_size (x)))
#else
	#define PART_MEM global
	#define MAX_SIZE(x) 
#endif


#define fast_mad(x, y, z) mad24(x, y, z)
#define fast_mul(x, y) mul24(x, y)








