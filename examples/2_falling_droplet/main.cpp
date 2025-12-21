// main.cpp
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include "particle.hpp"
#include "particle_filler.hpp"
#include "sph_solver.hpp"
#include "vector2d.hpp"
#include "vtk_io.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace fs = std::filesystem;

int main() {
    // Simulation parameters (change these without affecting logic)
    const int N = 2000;                     // number of particles
    const double density = 1000.0;          // reference density
    const double cs = 10.0;                 // sound speed
    const double viscosity = 0.005;         // dynamic viscosity
    const double sigma = 0.072;             // surface tension
    const Vector2D gravity = {0.0, -9.81};  // gravity acceleration

    const double box_length = 0.005;  // boundary box size

    const double radius = 0.001;                  // fill radius
    const Vector2D center = {0.0, 1.0 * radius};  // fill center

    const double dp = radius * std::sqrt(M_PI / N);  // particle spacing
    const double h = 1.3 * dp;                       // smoothing length

    const double CFL = 0.5;       // Courant number
    const double sim_time = 0.1;  // total simulation time in seconds
    const int frames = 1000;      // number of output frames

    const Vector2D initial_vel = {0.0, -0.5};  // initial velocity for all particles

    // Derived parameters
    const double max_dt = CFL * h / cs;
    const int min_num_steps = static_cast<int>(std::ceil(sim_time / max_dt));
    const int save_interval = static_cast<int>(std::ceil(static_cast<double>(min_num_steps) / frames));
    const int num_steps = frames * save_interval;
    const double dt = (num_steps == 0) ? 0.0 : sim_time / num_steps;

    // Create solver with parameters
    SPHSolver solver(h, dt, density, cs, viscosity, sigma, gravity, box_length);

    std::cout << "Kernel radius = " << h << std::endl;
    std::cout << "Time step = " << dt << std::endl;

    // Particles container
    std::vector<particle> particles;

    // Fill particles
    fill_circle(center, radius, particles, N, density, h);
    // fill_ellipse(center, radius * 1.5, radius / 1.5, particles, N, density, h);

    // Set initial velocity
    for (auto& p : particles) {
        p.vel = initial_vel;
    }

    // Output directory and initial save
    fs::create_directory("output");
    save_particles_to_vtk(particles, "output/particles_initial.vtk", 0);

    // Energy file
    // std::ofstream file_energy("energy.txt");

    // Timing
    auto start_time = std::chrono::high_resolution_clock::now();

    // Initial leapfrog setup
    solver.compute_forces(particles);
    for (auto& p : particles) {
        Vector2D a = p.force / p.mass;
        p.vel -= a * (dt / 2.0);
    }

    // Simulation loop
    for (int step = 0; step < num_steps; ++step) {
        solver.simulate_step(particles);

        // Save at intervals
        if (step % save_interval == 0) {
            std::ostringstream filename;
            filename << "output/particles_" << std::setw(6) << std::setfill('0') << (step / save_interval) << ".vtk";
            save_particles_to_vtk(particles, filename.str(), step);

            std::cout << "Time: " << step * dt << "/" << sim_time << std::endl;

            // Energy calculation
            /*
            double kineticEnergy = 0.0;
            double potentialEnergy = 0.0;
            for (const auto& p : particles) {
                kineticEnergy += 0.5 * p.mass * (p.vel.x * p.vel.x + p.vel.y * p.vel.y);
                potentialEnergy += p.mass * (-gravity.y * p.pos.y);
            }
            std::cout << "Kinetic Energy = " << kineticEnergy << std::endl;
            std::cout << "Potential Energy = " << potentialEnergy << std::endl;
            std::cout << "Total Energy = " << kineticEnergy + potentialEnergy << std::endl;
            file_energy << kineticEnergy << '\t' << potentialEnergy << '\t' << kineticEnergy + potentialEnergy << std::endl;
            */
        }
    }

    // file_energy.close();

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    std::cout << "Simulation complete. Total time: " << duration.count() / 1000.0 << " s" << std::endl;

    return 0;
}