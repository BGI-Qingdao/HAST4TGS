.PHONY: all clean

all: classify  classify_onelib

classify : classify.cpp gzstream/gzstream.C gzstream/gzstream.h
	g++ -c -g  gzstream/gzstream.C -I./gzstream -lz -o gzstream.o
	g++ -g -std=c++11 classify.cpp gzstream.o -lz -lpthread -o classify

classify_onelib : classify_onelib.cpp gzstream/gzstream.C gzstream/gzstream.h
	g++ -c -g  gzstream/gzstream.C -I./gzstream -lz -o gzstream.o
	g++ -g -std=c++11 classify_onelib.cpp gzstream.o -lz -lpthread -o classify_onelib

clean :
	rm -f classify
	rm -f *.o gzstream/*.o
