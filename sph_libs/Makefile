# Makefile for sph_libs (place this in mySPH/sph_libs)

CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -Werror -O3
INCLUDES = 
LIBS = 

OBJS = sph_kernel.o neighbor_finder.o force_calculator.o sph_solver.o particle_filler.o vtk_io.o

ifeq ($(OS),Windows_NT)
    RM = del /Q
else
    RM = rm -f
endif

all: $(OBJS)

sph_kernel.o: sph_kernel.cpp sph_kernel.hpp vector2d.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c sph_kernel.cpp

neighbor_finder.o: neighbor_finder.cpp neighbor_finder.hpp vector2d.hpp particle.hpp sph_kernel.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c neighbor_finder.cpp

force_calculator.o: force_calculator.cpp force_calculator.hpp vector2d.hpp particle.hpp sph_kernel.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c force_calculator.cpp

sph_solver.o: sph_solver.cpp sph_solver.hpp vector2d.hpp particle.hpp sph_kernel.hpp neighbor_finder.hpp force_calculator.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c sph_solver.cpp

particle_filler.o: particle_filler.cpp particle_filler.hpp vector2d.hpp particle.hpp sph_kernel.hpp neighbor_finder.hpp force_calculator.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c particle_filler.cpp

vtk_io.o: vtk_io.cpp vtk_io.hpp vector2d.hpp particle.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c vtk_io.cpp

clean:
	$(RM) $(OBJS)