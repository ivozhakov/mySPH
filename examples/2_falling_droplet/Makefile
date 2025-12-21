#Makefile for SPH simulation example(place this in mySPH / examples / <task_folder>)

CXX = g++ CXXFLAGS = -std = c++ 17 - Wall - Wextra - Werror - O3 LIB_DIR =../../ sph_libs INCLUDES = -I$(LIB_DIR) LIBS =

                                                                              OBJS = main.o LIB_OBJS = $(LIB_DIR) / sph_kernel.o $(LIB_DIR) / neighbor_finder.o $(LIB_DIR) /
                                                                                                       force_calculator.o $(LIB_DIR) / sph_solver.o $(LIB_DIR) / particle_filler.o $(LIB_DIR) /
                                                                                                       vtk_io.o

                                                                                                           EXEC = start.exe

                                                                                                                      ifeq($(OS), Windows_NT) RM = del / Q else RM = rm -
                                                                                                                                                                         f endif

                                                                                                                                                                             all : $(EXEC)

                                                                                                                                                                                       $(EXEC)
    : $(OBJS) $(LIB_OBJS) $(CXX) $(CXXFLAGS) $(INCLUDES) -
      o $(EXEC) $(OBJS) $(LIB_OBJS) $(LIBS)

              main.o : main.cpp $(LIB_DIR) /
          vector2d.hpp $(LIB_DIR) / particle.hpp $(LIB_DIR) / sph_kernel.hpp $(LIB_DIR) / neighbor_finder.hpp $(LIB_DIR) / force_calculator.hpp $(LIB_DIR) / sph_solver.hpp $(LIB_DIR) /
          particle_filler.hpp $(LIB_DIR) / vtk_io.hpp $(CXX) $(CXXFLAGS) $(INCLUDES) -
      c main.cpp

              $(LIB_DIR) /
          sph_kernel.o : $(LIB_DIR) / sph_kernel.cpp $(LIB_DIR) / sph_kernel.hpp $(LIB_DIR) / vector2d.hpp $(CXX) $(CXXFLAGS) $(INCLUDES) -
      c $ < -o $ @

                                                                                                                                                                             $(LIB_DIR) /
                                                                                                                                                                             neighbor_finder.o
    : $(LIB_DIR) / neighbor_finder.cpp $(LIB_DIR) / neighbor_finder.hpp $(LIB_DIR) / vector2d.hpp $(LIB_DIR) / particle.hpp $(LIB_DIR) / sph_kernel.hpp $(CXX) $(CXXFLAGS) $(INCLUDES) -
                                                                                                                                                                         c $ <
                                                                                                                                                                     -o $ @

                                                                                                                                                                             $(LIB_DIR) /
                                                                                                                                                                             force_calculator.o
    : $(LIB_DIR) / force_calculator.cpp $(LIB_DIR) / force_calculator.hpp $(LIB_DIR) / vector2d.hpp $(LIB_DIR) / particle.hpp $(LIB_DIR) / sph_kernel.hpp $(CXX) $(CXXFLAGS) $(INCLUDES) -
                                                                                                                                                                         c $ <
                                                                                                                                                                     -o $ @

                                                                                                                                                                             $(LIB_DIR) /
                                                                                                                                                                             sph_solver.o
    : $(LIB_DIR) / sph_solver.cpp $(LIB_DIR) / sph_solver.hpp $(LIB_DIR) / vector2d.hpp $(LIB_DIR) / particle.hpp $(LIB_DIR) / sph_kernel.hpp $(LIB_DIR) / neighbor_finder.hpp $(LIB_DIR) /
                                                                                                                                                                             force_calculator.hpp $(CXX)
                                                                                                                                                                                 $(CXXFLAGS)
                                                                                                                                                                                     $(INCLUDES) -
                                                                                                                                                                         c $ <
                                                                                                                                                                     -o $ @

                                                                                                                                                                             $(LIB_DIR) /
                                                                                                                                                                             particle_filler.o
    : $(LIB_DIR) / particle_filler.cpp $(LIB_DIR) / particle_filler.hpp $(LIB_DIR) / vector2d.hpp $(LIB_DIR) / particle.hpp $(LIB_DIR) / sph_kernel.hpp $(LIB_DIR) / neighbor_finder.hpp $(LIB_DIR) /
                                                                                                                                                                             force_calculator.hpp $(CXX)
                                                                                                                                                                                 $(CXXFLAGS)
                                                                                                                                                                                     $(INCLUDES) -
                                                                                                                                                                         c $ <
                                                                                                                                                                     -o $ @

                                                                                                                                                                             $(LIB_DIR) /
                                                                                                                                                                             vtk_io.o
    : $(LIB_DIR) / vtk_io.cpp $(LIB_DIR) / vtk_io.hpp $(LIB_DIR) / vector2d.hpp $(LIB_DIR) / particle.hpp $(CXX) $(CXXFLAGS) $(INCLUDES) - c $ < -o $ @

                                                                                                                                                                     clean : $(RM) $(OBJS) $(EXEC)