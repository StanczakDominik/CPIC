CC=g++
CFLAGS=-I.
DEPS = grid.hpp species.hpp
OBJ = grid.cpp species.hpp

all: main
main: main.o grid.o species.o
	$(CC) -o main main.o grid.o species.o $(CFLAGS) $(LIBS)

grid.o: grid.cpp grid.hpp species.hpp
	$(CC) -c grid.cpp $(CFLAGS)

species.o: species.cpp grid.hpp species.hpp
	$(CC) -c species.cpp $(CFLAGS)

clean:
	rm *.o
