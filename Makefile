CC=g++
CFLAGS=-I. -Wall -Werror -std=c++11 -lfftw3
CFILES=species.cpp grid.cpp simulation.cpp
HFILES=species.hpp grid.hpp simulation.hpp
OFILES=species.o grid.o simulation.o 

all: main
main: $(OFILES) main.o
	$(CC) -o main main.o $(OFILES) $(CFLAGS) 

main.o: main.cpp grid.o species.o simulation.o
	$(CC) -c main.cpp $(CFLAGS)
grid.o: grid.cpp species.hpp grid.hpp
	$(CC) -c grid.cpp $(CFLAGS)
species.o: species.cpp grid.hpp species.hpp
	$(CC) -c species.cpp $(CFLAGS)
simulation.o: $(CFILES) $(HFILES)
	$(CC) -c simulation.cpp $(CFLAGS)

run: all
	./main
clean:
	rm *.o
