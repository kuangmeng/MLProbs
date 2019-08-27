# profile -pg on compile & link
# %.o : %.cc ; $(CXX) -O3 -pg -D_FILE_OFFSET_BITS=64 -Wall -W  -c $<
%.o : %.cc ; $(CXX) -O3 -g -D_FILE_OFFSET_BITS=64 -Wall -W  -c $<
%.o : %.cpp ; $(CXX) -O3 -g -DNDEBUG -D_FILE_OFFSET_BITS=64 -Wall -W  -c $<

qscore:	\
	clineshift.o \
	comparemap.o \
	comparemsa.o \
	comparepair.o \
	fasta.o \
	fastq.o \
	gapscore.o \
	gapscore2.o \
	main.o \
	msa.o \
	perseq.o \
	options.o \
	qscore.o \
	seq.o \
	sab.o \
	sumpairs.o \
	tc.o \
	usage.o \
	utils.o \

# profile needs -pg
#	$(CXX)  -O3 -pg -o qscore $^
	$(CXX)  -O3 -o qscore $^
