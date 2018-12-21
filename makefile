CC = g++
CXXFLAGS = -g -O3 -lm -lboost_program_options -std=c++11 

BayesC: 
	$(CC) -o BayesC BayesC.cpp BayesC_distributions.cpp Sampling_functions.cpp $(CXXFLAGS)

clean: 
	rm BayesC