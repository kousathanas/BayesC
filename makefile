CC = g++
CXXFLAGS = -g -O3 -lm -lgsl -std=c++11

BayesC: 
	$(CC) -o BayesC BayesC.cpp  $(CXXFLAGS)

clean: 
	rm BayesC