// particle_filler.cpp
#include "particle_filler.hpp"

#include <cmath>
#include <iostream>
#include <random>

#include "force_calculator.hpp"
#include "neighbor_finder.hpp"
#include "sph_kernel.hpp"

void fill_circle(Vector2D center, double radius, std::vector<particle>& particles, int N, double density0, double h) {
    if (N <= 0 || radius <= 0.0) return;

    particles.clear();
    particles.reserve(N);

    double area = M_PI * radius * radius;
    double area_per_particle = area / N;
    double spacing = std::sqrt(area_per_particle) * 1.071;

    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-0.2 * spacing, 0.2 * spacing);

    double hex_offset = spacing * std::sqrt(3.0) / 2.0;
    int num_rows = static_cast<int>((2.0 * radius + 2.0 * spacing) / hex_offset) + 2;

    int particles_added = 0;
    double buffer = spacing;
    for (int row = 0; row < num_rows && particles_added < N; ++row) {
        double y = center.y - radius - buffer + row * hex_offset;
        double x_offset = (row % 2 == 0) ? 0.0 : 0.5 * spacing;

        int num_cols = static_cast<int>((2.0 * radius + 2.0 * buffer) / spacing) + 2;
        for (int col = 0; col < num_cols && particles_added < N; ++col) {
            double x = center.x - radius - buffer + col * spacing + x_offset;

            double dx = x - center.x;
            double dy = y - center.y;
            if (dx * dx + dy * dy <= radius * radius) {
                double jitter_x = dis(gen);
                double jitter_y = dis(gen);
                particles.emplace_back(x + jitter_x, y + jitter_y);
                particles_added++;
            }
        }
    }

    if (particles_added < N) {
        std::cout << "Warning: Added only " << particles_added << " particles out of " << N << std::endl;
    }

    find_neighbors(particles, h);
    compute_density(particles, density0);

    double avg_rho = 0.0;
    for (const auto& p : particles) {
        avg_rho += p.density;
    }
    avg_rho /= particles.size();

    double scale = density0 / avg_rho;
    for (auto& p : particles) {
        p.mass *= scale;
    }
    compute_density(particles, density0);
}

void fill_ellipse(Vector2D center, double rx, double ry, std::vector<particle>& particles, int N, double density0, double h) {
    if (N <= 0 || rx <= 0.0 || ry <= 0.0) return;

    particles.clear();
    particles.reserve(N);

    double area = M_PI * rx * ry;
    double area_per_particle = area / N;
    double spacing = std::sqrt(area_per_particle) * 1.071;

    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-0.2 * spacing, 0.2 * spacing);

    double hex_offset = spacing * std::sqrt(3.0) / 2.0;
    int num_rows = static_cast<int>((2.0 * ry + 2.0 * spacing) / hex_offset) + 2;

    int particles_added = 0;
    double buffer = spacing;
    for (int row = 0; row < num_rows && particles_added < N; ++row) {
        double y = center.y - ry - buffer + row * hex_offset;
        double x_offset = (row % 2 == 0) ? 0.0 : 0.5 * spacing;

        int num_cols = static_cast<int>((2.0 * rx + 2.0 * buffer) / spacing) + 2;
        for (int col = 0; col < num_cols && particles_added < N; ++col) {
            double x = center.x - rx - buffer + col * spacing + x_offset;

            double dx = x - center.x;
            double dy = y - center.y;
            if ((dx * dx) / (rx * rx) + (dy * dy) / (ry * ry) <= 1.0) {
                double jitter_x = dis(gen);
                double jitter_y = dis(gen);
                particles.emplace_back(x + jitter_x, y + jitter_y);
                particles_added++;
            }
        }
    }

    if (particles_added < N) {
        std::cout << "Warning: Added only " << particles_added << " particles out of " << N << std::endl;
    }

    find_neighbors(particles, h);
    compute_density(particles, density0);

    double avg_rho = 0.0;
    for (const auto& p : particles) {
        avg_rho += p.density;
    }
    avg_rho /= particles.size();

    double scale = density0 / avg_rho;
    for (auto& p : particles) {
        p.mass *= scale;
    }
    compute_density(particles, density0);
}