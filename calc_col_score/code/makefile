CC     = g++
CFLAGS = -c -Wall -O3
LFLAGS = -g -Wall -O3 

all: calc_col_score

install: calc_col_score

clean : 
	rm -f *.o *.exe

OBJECTS = msa_functions.o count_functions.o set_functions.o 


calc_col_score : $(OBJECTS) defs.h calc_col_score.o
	$(CC) $(LFLAGS) -o $@ $(OBJECTS) $@.o 


%.o : %.cpp defs.h
	$(CC) $(CFLAGS) $<
                
