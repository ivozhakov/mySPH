// particle_filler.hpp
#ifndef PARTICLE_FILLER_HPP
#define PARTICLE_FILLER_HPP

#include <vector>

#include "particle.hpp"
#include "vector2d.hpp"

void fill_circle(Vector2D center, double radius, std::vector<particle>& particles, int N, double density0, double h);
void fill_ellipse(Vector2D center, double rx, double ry, std::vector<particle>& particles, int N, double density0, double h);

#endif  // PARTICLE_FILLER_HPP