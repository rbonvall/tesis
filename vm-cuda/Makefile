NVCC = nvcc
CC = g++
CFLAGS = -Wall -I ~/codigo/boost_1_42_0
TARGETS = init vm particles

all: $(TARGETS)

init.o: init.cpp vm.hpp
	$(CC) $(CFLAGS) -c $<
init: init.o
	$(CC) $< -o $@

vm.o: vm.cpp vm.hpp
	$(CC) $(CFLAGS) -c $< 
vm: vm.o
	$(CC) $< -o $@

particles: init
	./$< > $@

clean:
	rm -f *.o $(TARGETS)

.PHONY: all clean
