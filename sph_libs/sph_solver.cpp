// sph_solver.cpp
#include "sph_solver.hpp"

#include "force_calculator.hpp"
#include "neighbor_finder.hpp"
#include "sph_kernel.hpp"

SPHSolver::SPHSolver(double smoothing_length, double time_step, double ref_density, double pressure_stiffness, double fluid_viscosity, double surface_tension, Vector2D acceleration, double box_size)
    : h(smoothing_length), dt(time_step), density0(ref_density), stiffness(pressure_stiffness), viscosity(fluid_viscosity), sigma(surface_tension), gravity(acceleration), box_length(box_size) {
    SPHKernel::initialize(h);
}

void SPHSolver::compute_forces(std::vector<particle>& particles) {
    for (auto& p : particles) {
        p.force = {0.0, 0.0};
    }
    find_neighbors(particles, h);
    compute_pressure_forces(particles, stiffness, density0, h);
    compute_viscous_forces(particles, viscosity, h);
    compute_surface_tension_CSS(particles, sigma, h);
    for (auto& p : particles) {
        p.force += gravity * p.mass;
    }
}

void SPHSolver::simulate_step(std::vector<particle>& particles) {
    for (auto& p : particles) {
        p.vel += (p.force / p.mass) * dt;
        p.pos += p.vel * dt;
    }
    apply_boundaries(particles);
    compute_forces(particles);
}

double SPHSolver::get_smoothing_length() const { return h; }
double SPHSolver::get_time_step() const { return dt; }

void SPHSolver::apply_boundaries(std::vector<particle>& particles) const {
    for (auto& p : particles) {
        double k = 1.0;
        if (p.pos.y < 0) {
            double delta = 0 - p.pos.y;
            p.pos.y = 0 + delta;
            p.vel.y = -p.vel.y * k;
        }

        if (p.pos.y > box_length) {
            double delta = p.pos.y - box_length;
            p.pos.y = box_length - delta;
            p.vel.y = -p.vel.y * k;
        }

        if (p.pos.x > box_length) {
            double delta = p.pos.x - box_length;
            p.pos.x = box_length - delta;
            p.vel.x = -p.vel.x * k;
        }

        if (p.pos.x < -box_length) {
            double delta = -box_length - p.pos.x;
            p.pos.x = -box_length + delta;
            p.vel.x = -p.vel.x * k;
        }
    }
}