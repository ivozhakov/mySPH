// force_calculator.hpp
#ifndef FORCE_CALCULATOR_HPP
#define FORCE_CALCULATOR_HPP

#include <vector>

#include "particle.hpp"
#include "vector.hpp"

template<int Dim>
void compute_surface_tension_CSS(std::vector<Particle<Dim>>& particles, double sigma, double h);

template<int Dim>
void compute_density(std::vector<Particle<Dim>>& particles, double density0);

template<int Dim>
void density_normalization(std::vector<Particle<Dim>>& particles, double density0);

template<int Dim>
void compute_pressure_forces(std::vector<Particle<Dim>>& particles, double cs, double density0, double h);

template<int Dim>
void compute_viscous_forces(std::vector<Particle<Dim>>& particles, double viscosity, double h);

#endif  // FORCE_CALCULATOR_HPP