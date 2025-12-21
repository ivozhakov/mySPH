// neighbor_finder.hpp
#ifndef NEIGHBOR_FINDER_HPP
#define NEIGHBOR_FINDER_HPP

#include <vector>

#include "particle.hpp"
#include "vector2d.hpp"

void find_neighbors(std::vector<particle>& particles, double h);

#endif  // NEIGHBOR_FINDER_HPP