// sph_solver.cpp
#include "sph_solver.hpp"

#include "force_calculator.hpp"
#include "neighbor_finder.hpp"
#include "sph_kernel.hpp"

template<int Dim>
SPHSolver<Dim>::SPHSolver(double smoothing_length
                        , double time_step
                        , double ref_density
                        , double sound_velocity
                        , double fluid_viscosity
                        , double surface_tension
                        , Vector<Dim> acceleration
                        , std::pair<Vector<Dim>, Vector<Dim>> box)
    : h(smoothing_length)
    , dt(time_step)
    , density0(ref_density)
    , cs(sound_velocity)
    , viscosity(fluid_viscosity)
    , sigma(surface_tension)
    , gravity(acceleration)
    , box(box) {
    SPHKernel::initialize<Dim>(h);
}

template<int Dim>
void SPHSolver<Dim>::compute_forces(std::vector<Particle<Dim>>& particles) {
    for (auto& p : particles) {
        p.force = Vector<Dim>();
    }
    find_neighbors(particles, h);
    compute_pressure_forces(particles, cs, density0, h);
    compute_viscous_forces(particles, viscosity, h);
    compute_surface_tension_CSS(particles, sigma, h);
    for (auto& p : particles) {
        p.force += gravity * p.mass;
    }
}

template<int Dim>
void SPHSolver<Dim>::simulate_step(std::vector<Particle<Dim>>& particles) {
    for (auto& p : particles) {
        p.vel += (p.force / p.mass) * dt;
        p.pos += p.vel * dt;
    }
    apply_boundaries(particles);
    compute_forces(particles);
}

template<int Dim>
double SPHSolver<Dim>::get_smoothing_length() const { return h; }

template<int Dim>
double SPHSolver<Dim>::get_time_step() const { return dt; }

template<int Dim>
void SPHSolver<Dim>::apply_boundaries(std::vector<Particle<Dim>>& particles) const {
    for (auto& p : particles) {
        double k = 1.0;  // restitution coefficient
        for (int d = 0; d < Dim; ++d) {
            double min_d = box.first[d];
            double max_d = box.second[d];
            if (p.pos[d] < min_d) {
                double delta = min_d - p.pos[d];
                p.pos[d] = min_d + delta;
                p.vel[d] = -p.vel[d] * k;
            }
            if (p.pos[d] > max_d) {
                double delta = p.pos[d] - max_d;
                p.pos[d] = max_d - delta;
                p.vel[d] = -p.vel[d] * k;
            }
        }
    }
}

template class SPHSolver<2>;
template class SPHSolver<3>;