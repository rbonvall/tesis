NVCC = nvcc
CC = g++
CXXFLAGS = -Wall $(INCLUDEFLAGS)

ifeq ($(shell hostname),gp01)
    BOOST = ~/codigo/boost_1_42_0
    BOOSTHEADERS = $(BOOST)
    BOOSTLIB = $(BOOST)/stage/lib/

    CUDA = /opt/cuda
    CUDAHEADERS = $(CUDA)/include
else
    CUDA = /usr/local/cuda
    CUDAHEADERS = $(CUDA)/include
    NVCCFLAGS = #-deviceemu
    BOOSTHEADERS = /usr/include/boost
    LDFLAGS += -L$(CUDA)/lib64 -L$(CUDA)/lib
    #LD_LIBRARY_PATH := /usr/local/cuda/lib64:/usr/local/cuda/lib:$(LD_LIBRARY_PATH)
    #LIBRARY_PATH = $(LD_LIBRARY_PATH)
endif

THRUSTHEADERS = thrust

INCLUDEFLAGS = -I $(BOOSTHEADERS) -I $(CUDAHEADERS) -I $(THRUSTHEADERS)
LDFLAGS += -lm -lboost_program_options -lgomp -lcudart

SRC = $(wildcard *.cpp)
TARGETS = init analytical-vel vel-eval integrate integrate-gpu

all: $(TARGETS)
init: init.o vm.o
analytical-vel: analytical-vel.o vm.o
vel-eval: vel-eval.o vm.o
integrate: integrate.o vm.o
integrate-gpu: cuda-integrate.o vm.o

cuda-integrate.o: cuda-integrate.cu cuda-integrate.hpp
	$(NVCC) $(NVCCFLAGS) $(INCLUDEFLAGS) $< -c -o $@

%.d: %.cpp
	@set -e; rm -f $@; \
	$(CC) -M $(CXXFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

-include $(SRC:.cpp=.d)

clean:
	rm -f *.o $(TARGETS) *.d

.PHONY: all clean
