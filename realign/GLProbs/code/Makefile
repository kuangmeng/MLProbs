
CXXOBJS = MSA.o MSAGuideTree.o MSAClusterTree.o MSAPartProbs.o MSAReadMatrix.o main.o

OPENMP = -fopenmp
CXX = g++
COMMON_FLAGS = -O3 -lm $(OPENMP) -Wall -funroll-loops -I . -I /usr/include
CXXFLAGS = $(COMMON_FLAGS)

EXEC = glprobs

all: $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -o $(EXEC) $(CXXOBJS) $(NVCCOBJS) $(NVCCLIBS)
	strip $(EXEC)
clean:
	rm -rf *.o $(EXEC)

