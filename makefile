CC = mpic++
CXXFLAGS = -g -O3 -lm -lboost_program_options -std=c++11 -fopenmp -lmpi -I/usr/lib/openmpi/include 

BayesC: 
	$(CC) -o BayesC BayesC.cpp data.cpp mympi.cpp gadgets.cpp $(CXXFLAGS)

clean: 
	rm BayesC