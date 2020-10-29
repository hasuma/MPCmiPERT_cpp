-include subdir.mk

CXX := g++
CXX_FLAGS := -Wall -Wextra -std=c++17 -ggdb
CPPFLAGS = -Wall -pedantic -Wshadow -Wfloat-equal -O3 -Wconversion -Wsign-conversion -fPIC -DLINUX -D__USE_LONG_INTEGERS__ -D__USE_LONG_FINTS__ -D__NO_COPYRIGHT__ 

BIN := bin
SRC := src
INCLUDE :=

LIBRARIES :=-larmadillo -llapack -lblas -L/include/qpOASES/bin -lqpOASES
EXECUTABLE := bin/executavel
QPOASES_OBJECTS = \
	include/qpOASES/src/SQProblem.o \
	include/qpOASES/src/QProblem.o \
	include/qpOASES/src/QProblemB.o \
	include/qpOASES/src/SQProblemSchur.o\
	include/qpOASES/src/Bounds.o \
	include/qpOASES/src/Constraints.o \
	include/qpOASES/src/SubjectTo.o\
	include/qpOASES/src/Indexlist.o \
	include/qpOASES/src/Flipper.o\
	include/qpOASES/src/Utils.o\
	include/qpOASES/src/Options.o\
	include/qpOASES/src/Matrices.o \
	include/qpOASES/src/MessageHandling.o \
	include/qpOASES/src/SparseSolver.o

all:	$(EXECUTABLE)

run:	
	clear	
	@echo "Executing..."
	./$(EXECUTABLE)

$(EXECUTABLE): 	$(OBJS)
		@echo "Building..."
		$(CXX) $(CPPFLAGS) -o $@ $(QPOASES_OBJECTS) "src/$<" $(LIBRARIES)
		@echo 'Finished building target: $@'
clean:
	@echo "Clearing..."
	rm -rf src/$(OBJS)