CXX=g++
STD=-std=c++17
WAR=-Wall -Wextra -pedantic
BUG=-g
OPT=
CXX_FLAGS=$(BUG) $(WAR) $(OPT) $(STD)
.PHONY: all clean

TARGETS=NSim
OBJECTS=Source.o NSim.o Integrator.o PoissonSolver.o MassTree.o Particle.o

# A note on variables:
# $@: the target filename.
# $*: the target filename without the file extension.
# $<: the first prerequisite filename.
# $^: the filenames of all the prerequisites, separated by spaces, discard duplicates.
# $+: similar to $^, but includes duplicates.
# $?: the names of all prerequisites that are newer than the target, separated by spaces.

all: $(TARGETS)

%.o: %.cpp
	$(CXX) $(CXX_FLAGS) -c $<

Source.o: Source.cpp
NSim.o: NSim.cpp NSim.h
Integrator.o: Integrator.cpp Integrator.h
PoissonSolver.o: PoissonSolver.cpp PoissonSolver.h
MassTree.o: MassTree.cpp MassTree.h
Particle.o: Particle.cpp Particle.h

NSim: $(OBJECTS)
	$(CXX) $(CXX_FLAGS) -o $@ $^
#	rm -f *.o

clean:
	rm -f *.o *~ $(TARGETS)
