NVCC = nvcc
CC = g++
CFLAGS = -Wall -I ~/codigo/boost_1_42_0
TARGETS = init particles

all: $(TARGETS)

init.o: init.cpp vm.o
	$(CC) $(CFLAGS) -c $<
init: init.o
	$(CC) $< -o $@

vm.o: vm.cpp vm.hpp
	$(CC) $(CFLAGS) -c $< 

particles: init
	./$< > $@

clean:
	rm -f *.o $(TARGETS)

.PHONY: all clean
