NVCC = nvcc
CC = g++
INCLUDEFLAGS = -I ~/codigo/boost_1_42_0
CXXFLAGS = -Wall $(INCLUDEFLAGS)
LDFLAGS += -lm -lboost_program_options
TARGETS = init particles vel-eval

SRC = init.cpp vel-eval.cpp vm.cpp
#MAKEDEPEND = g++ -M $(CPPFLAGS) -o $*.d $<

all: $(TARGETS)
vel-eval: vel-eval.o vm.o
init: init.o vm.o

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
