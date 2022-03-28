CC=gcc
CXX=g++
RM=rm -f
CPPFLAGS=-g -std=c++11
LDFLAGS=-g -std=c++11
LDLIBS=-g -std=c++11 # -O3

SRCS=EulerDiffusionSolver.cpp DiffusionProblem.cpp AbstractVector.cpp Vector.cpp SparseVector.cpp SparseMatrix.cpp globals.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

runtests: test
	./test

run: main
	./main

all: test main

test: test.o $(OBJS)
	$(CXX) $(LDFLAGS) -o test test.o $(OBJS) $(LDLIBS)

main: main.o $(OBJS)
	$(CXX) $(LDFLAGS) -o main main.o $(OBJS) $(LDLIBS)

main.o: main.cpp $(SRCS)

test.o: test.cpp $(SRCS)


SparseMatrix.o: SparseMatrix.cpp SparseMatrix.hpp AbstractVector.hpp SparseVector.hpp Vector.hpp

Vector.o: Vector.cpp Vector.hpp

clean:
	$(RM) $(OBJS) *.o *.$(DLEXT)
	$(RM) test
