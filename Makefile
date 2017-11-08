CC=gcc
CFLAGS=-I.
DEPS = grid.hpp species.hpp
OBJ = grid.cpp species.hpp

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

