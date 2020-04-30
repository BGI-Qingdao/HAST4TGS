.PHONY: all clean

all: classify 

classify : classify.cpp gzstream/gzstream.C gzstream/gzstream.h
	g++ -c -g  gzstream/gzstream.C -I./gzstream -lz -o gzstream.o
	g++ -g -std=c++11 classify.cpp gzstream.o -lz -lpthread -o classify

clean :
	rm -f classify
	rm -f *.o gzstream/*.o
