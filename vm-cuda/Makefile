NVCC = nvcc
CC = g++
CXXFLAGS = -Wall $(INCLUDEFLAGS)

BOOST = ~/codigo/boost_1_42_0
BOOSTHEADERS = $(BOOST)
BOOSTLIB = $(BOOST)/stage/lib/

INCLUDEFLAGS = -I $(BOOSTHEADERS)
LDFLAGS += -lm -lboost_program_options -lgomp

SRC = $(wildcard *.cpp)
TARGETS = init analytical-vel vel-eval particles

all: $(TARGETS)
init: init.o vm.o
analytical-vel: analytical-vel.o vm.o
vel-eval: vel-eval.o vm.o

%.d: %.cpp
	@set -e; rm -f $@; \
	$(CC) -M $(CXXFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

-include $(SRC:.cpp=.d)

particles: init
	./$< > $@
clean:
	rm -f *.o $(TARGETS) *.d

.PHONY: all clean
