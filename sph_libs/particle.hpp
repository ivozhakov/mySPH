// particle.hpp
#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <vector>

#include "vector.hpp"

template<int Dim>
struct Particle {
    Vector<Dim> pos;
    Vector<Dim> vel;
    Vector<Dim> force;
    double density;
    double pressure;
    double mass;
    double c;
    Vector<Dim> grad_c;
    double mag_grad_c;
    double kappa;
    double cspm_den;  // For CSPM normalization

    std::vector<size_t> neighbors;  // Neighbor lists
    // std::vector<size_t> neighbors_CSS;  // For CSS if needed

    Particle() : density(1000.0), pressure(0.0), mass(0.02), c(0.0), mag_grad_c(0.0), kappa(0.0), cspm_den(0.0) {}
};

#endif  // PARTICLE_HPP