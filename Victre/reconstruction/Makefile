################################################################################
#
# This program compiles the C codes for DBT reconstruction.
# Usage: make
 
################################################################################

################################################################################
#
# Build script for project
#
################################################################################

# Add source files here
EXECUTABLE	:= ./FBP
# C/C++ source files (compiled with gcc / c++)
#CCFILES		:= oclDeviceQuery.cpp


################################################################################
# Rules and targets

# Compilers
CXX        := g++
LINK       := g++ -fPIC
CFLAGS     := -O3 -w

default: conebeam
	$(CXX) -o $(EXECUTABLE) FBP.o -lm -lfftw3f -L lib

conebeam: FBP_DBTrecon.c
	$(CXX) $(CFLAGS) -o FBP.o -c FBP_DBTrecon.c

clean: 
	rm *.o $(EXECUTABLE);

