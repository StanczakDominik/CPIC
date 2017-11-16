CC=g++
CFLAGS=-I. -Wall -Werror -std=c++11 -lfftw3
CFILES=species.cpp grid.cpp simulation.cpp
HFILES=species.hpp grid.hpp simulation.hpp
OFILES=species.o grid.o simulation.o

all: main
main: main.cpp $(OFILES)
	$(CC) -o main main.o $(OFILES) $(CFLAGS) 

grid.o: $(CFILES) $(HFILES)
	$(CC) -c grid.cpp $(CFLAGS)
species.o: $(CFILES) $(HFILES)
	$(CC) -c species.cpp $(CFLAGS)
simulation.o: $(CFILES) $(HFILES)
	$(CC) -c simulation.cpp $(CFLAGS)

clean:
	rm *.o
