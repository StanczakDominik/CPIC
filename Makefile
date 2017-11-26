CC=g++
CFLAGS=-I. -Wall -Werror -Wextra -std=c++11 -lfftw3 -O2 -march=native
# HDF5_FLAGS= -D_FORTIFY_SOURCE=2 -march=x86-64 -mtune=generic -O2 -pipe -fstack-protector-strong -L/usr/lib -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 -Wl,-O1,--sort-common,--as-needed,-z,relro -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib
# CFLAGS += $(HDF5_FLAGS)

all: main
main: main.o grid.o species.o simulation.o temporal.o
	$(CC) -o main main.o grid.o species.o simulation.o temporal.o $(CFLAGS) 

main.o: main.cpp grid.o species.o simulation.o temporal.o
	$(CC) -c main.cpp $(CFLAGS)
simulation.o: simulation.cpp species.o grid.o temporal.o
	$(CC) -c simulation.cpp $(CFLAGS)
species.o: species.cpp temporal.o grid.o
	$(CC) -c species.cpp $(CFLAGS)
grid.o: grid.cpp temporal.o
	$(CC) -c grid.cpp $(CFLAGS)
temporal.o: temporal.cpp
	$(CC) -c temporal.cpp $(CFLAGS)

run: all
	./main
clean:
	rm *.o
