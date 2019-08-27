cmake_minimum_required (VERSION 2.6)
project (Alignment)

set(PAIRWISE_FILES
	Pairwise/Backtrack.cpp
	Pairwise/Backtrack.h
	Pairwise/AminoAcidMatrices.cpp
	Pairwise/AminoAcidMatrices.hpp
	Pairwise/AminoAcidMatrices.h
	Pairwise/IPairwiseAligner.h
	Pairwise/NeedlemanWunsch.cpp
	Pairwise/NeedlemanWunsch.h
	Pairwise/NucleotideMatrices.cpp
	Pairwise/NucleotideMatrices.h
	Pairwise/NucleotideMatrices.hpp
	Pairwise/NucleotideHmmGTR5.cpp
	Pairwise/NucleotideHmmGTR5.h
	Pairwise/PairHmm.h
	Pairwise/PairHmm.cpp
	Pairwise/PairHmmFactory.h
	Pairwise/PfamHmm5.cpp
	Pairwise/PfamHmm5.h
	Pairwise/ProteinHmm3.h
	Pairwise/ProteinHmm3.cpp
	Pairwise/ProteinHmm5.h
	Pairwise/ProteinHmm5.cpp
	Pairwise/Scoring.h
	Pairwise/ScoringFactory.h	
	Pairwise/SmithWaterman.cpp
	Pairwise/SmithWaterman.h)

set(DATA_STRUCTURES_FILES
	DataStructures/ContiguousMultiSequence.cpp
	DataStructures/ContiguousMultiSequence.h
	DataStructures/FilteredSparseMatrix.cpp
	DataStructures/FilteredSparseMatrix.h
	DataStructures/FilteredSparseMatrix.hpp
	DataStructures/ISequenceSet.h
	DataStructures/ISparseMatrix.h
	DataStructures/JaggedMatrix.h
	DataStructures/MultiSequence.h
	DataStructures/MultiSequence.cpp
	DataStructures/PackedSparseMatrix.cpp
	DataStructures/PackedSparseMatrix.h
	DataStructures/PackedSparseMatrix.hpp
	DataStructures/Sequence.cpp
	DataStructures/Sequence.h
	DataStructures/SequenceIO.h
	DataStructures/SequenceIO.cpp
	DataStructures/SparseEntry.h
	DataStructures/SparseHelper.cpp
	DataStructures/SparseHelper.h
	DataStructures/SparseMatrixType.h)	
	
set(MSA_POSTERIOR_FILES
	Multiple/BufferSet.h
	Multiple/IAlgorithmStage.h
	Multiple/IPartitionFunctionParams.h
	Multiple/ParallelProbabilisticModel.cpp
	Multiple/ParallelProbabilisticModel.h
	Multiple/ExpPartitionFunctionParams.h
	Multiple/LogPartitionFunctionParams.h
	Multiple/PartitionFunction.cpp
	Multiple/PartitionFunction.h
	Multiple/PosteriorStage.cpp
	Multiple/PosteriorStage.h
	Multiple/ProbabilisticParams.cpp
	Multiple/ProbabilisticParams.h
	Multiple/ProbabilisticModel.h
	Multiple/ProbabilisticModel.cpp

)

set(MSA_CONSISTENCY_FILES
	Multiple/ConsistencyStage.cpp
	Multiple/ConsistencyStage.h
	Multiple/RelaxationTask.cpp
	Multiple/Selectivity.h
)

set(MSA_FINAL_FILES
	Multiple/AlignmentEvaluatorBase.h
	Multiple/AminoAcidProperties.h
	Multiple/AminoAcidProperties.cpp
	Multiple/EntropyEvaluator.h
	Multiple/EntropyEvaluator.cpp
	Multiple/ColumnRefinement.cpp
	Multiple/ColumnRefinement.h
	Multiple/ConstructionStage.cpp
	Multiple/ConstructionStage.h
	Multiple/RefinementBase.cpp
	Multiple/RefinementBase.h
	Multiple/RandomRefinement.cpp
	Multiple/RandomRefinement.h
	Multiple/ScoringRefinement.cpp
	Multiple/ScoringRefinement.h
	Multiple/TreeRefinement.cpp
	Multiple/TreeRefinement.h
)

set(MSA_PHYLOGENY_FILES
	Multiple/Node.h
	Multiple/ClusterTree.cpp
	Multiple/ClusterTree.h
	Multiple/GuideTree.cpp
	Multiple/GuideTree.h
	Multiple/NewickTree.h
	Multiple/NewickTree.cpp
	Multiple/PhylipTree.h
	Multiple/PhylipTree.cpp
	Multiple/TreeGrammar.h
	Multiple/SLinkTree.cpp
	Multiple/SLinkTree.h
)

set(MSA_MAIN_FILES
	Multiple/AuxiliaryStructures.h	
	Multiple/BasicMSA.h
	Multiple/BasicMSA.cpp
	Multiple/Common.h	
	Multiple/Configuration.cpp
	Multiple/Configuration.h
	Multiple/ConfigurationHelper.h
	Multiple/ExtendedMSA.cpp
	Multiple/ExtendedMSA.h
	Multiple/ScoreType.h	
	)


include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})	
include_directories(${Extra_INCLUDE_DIRS})
	
#	
source_group("Pairwise" FILES ${PAIRWISE_FILES})
source_group("DataStructures" FILES ${DATA_STRUCTURES_FILES})
source_group("Multiple\\Consistency" FILES ${MSA_CONSISTENCY_FILES})
source_group("Multiple\\Final" FILES ${MSA_FINAL_FILES})
source_group("Multiple\\Phylogeny" FILES ${MSA_PHYLOGENY_FILES})
source_group("Multiple\\Posterior" FILES ${MSA_POSTERIOR_FILES})
source_group("Multiple" FILES ${MSA_MAIN_FILES})
source_group("" FILES Alignment.h)

set(ALL_FILES Alignment.h ${PAIRWISE_FILES} ${DATA_STRUCTURES_FILES} ${MSA_CONSISTENCY_FILES} ${MSA_FINAL_FILES} ${MSA_PHYLOGENY_FILES} ${MSA_POSTERIOR_FILES} ${MSA_MAIN_FILES})

add_library(Alignment STATIC ${ALL_FILES} )
target_link_libraries(Alignment 
	debug ${TCMalloc_LIBRARIES_DEBUG} 
	optimized ${TCMalloc_LIBRARIES_RELEASE}
	general Common Hardware ${System_LIBRARIES})
