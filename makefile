CC = g++
CXXFLAGS = -g -O3 -lm -lgsl -lboost_program_options -std=c++11

BayesC: 
	$(CC) -o BayesC BayesC.cpp  $(CXXFLAGS)

clean: 
	rm BayesC