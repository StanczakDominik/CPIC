CC=g++
CFLAGS=-I. -Wall -Werror -std=c++11 -lfftw3

all: main
main: main.o grid.o species.o simulation.o temporal.o
	$(CC) -o main main.o grid.o species.o simulation.o temporal.o $(CFLAGS) 

main.o: main.cpp grid.o species.o simulation.o temporal.o
	$(CC) -c main.cpp $(CFLAGS)
grid.o: temporal.o
	$(CC) -c grid.cpp $(CFLAGS)
species.o: temporal.o
	$(CC) -c species.cpp $(CFLAGS)
simulation.o: species.o grid.o temporal.o
	$(CC) -c simulation.cpp $(CFLAGS)
temporal.o: 
	$(CC) -c temporal.cpp $(CFLAGS)

run: all
	./main
clean:
	rm *.o
