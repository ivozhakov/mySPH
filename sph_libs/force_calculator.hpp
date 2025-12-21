// force_calculator.hpp
#ifndef FORCE_CALCULATOR_HPP
#define FORCE_CALCULATOR_HPP

#include <vector>

#include "particle.hpp"
#include "vector2d.hpp"

void compute_surface_tension_CSS(std::vector<particle>& particles, double sigma, double h);
void compute_density(std::vector<particle>& particles, double density0);
void density_normalization(std::vector<particle>& particles, double density0);
void compute_pressure_forces(std::vector<particle>& particles, double cs, double density0, double h);
void compute_viscous_forces(std::vector<particle>& particles, double viscosity, double h);

#endif  // FORCE_CALCULATOR_HPP