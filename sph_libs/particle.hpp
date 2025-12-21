// particle.hpp
#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <vector>

#include "vector2d.hpp"

struct particle {
    Vector2D pos;
    Vector2D vel;
    Vector2D force;
    double density;
    double pressure;
    double mass;
    double c;
    Vector2D grad_c;
    double mag_grad_c;
    double kappa;
    double cspm_den;  // For CSPM normalization

    std::vector<size_t> neighbors;  // Neighbor lists
    // std::vector<size_t> neighbors_CSS;  // For CSS if needed

    particle(double x = 0.0, double y = 0.0)
        : pos(x, y), vel(0.0, 0.0), force(0.0, 0.0), density(1000.0), pressure(0.0), mass(0.02), c(0.0), grad_c(0.0, 0.0), mag_grad_c(0.0), kappa(0.0), cspm_den(0.0) {}
};

#endif  // PARTICLE_HPP