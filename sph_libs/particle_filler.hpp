// particle_filler.hpp
#ifndef PARTICLE_FILLER_HPP
#define PARTICLE_FILLER_HPP

#include <vector>

#include "particle.hpp"
#include "vector.hpp"

template <int Dim>
void fill_ellipsoid(Vector<Dim> center, std::array<double, Dim> radii, std::vector<Particle<Dim>>& particles, int N, double density0, double h);

#endif  // PARTICLE_FILLER_HPP