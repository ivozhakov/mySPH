// neighbor_finder.hpp
#ifndef NEIGHBOR_FINDER_HPP
#define NEIGHBOR_FINDER_HPP

#include <vector>

#include "particle.hpp"
#include "vector.hpp"

template<int Dim>
void find_neighbors(std::vector<Particle<Dim>>& particles, double h);

#endif  // NEIGHBOR_FINDER_HPP