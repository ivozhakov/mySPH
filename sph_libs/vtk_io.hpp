// vtk_io.hpp
#ifndef VTK_IO_HPP
#define VTK_IO_HPP

#include <string>
#include <vector>

#include "particle.hpp"

void save_particles_to_vtk(const std::vector<particle>& particles, const std::string& filename, int timestep = 0);

#endif  // VTK_IO_HPP