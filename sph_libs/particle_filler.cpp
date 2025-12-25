// particle_filler.cpp
#include "particle_filler.hpp"

#include <cmath>
#include <iostream>
#include <random>

#include "force_calculator.hpp"
#include "neighbor_finder.hpp"
#include "sph_kernel.hpp"

template <int Dim>
void fill_ellipsoid(Vector<Dim> center, std::array<double, Dim> radii, std::vector<Particle<Dim>>& particles, int N, double density0, double h) {
    particles.clear();
    particles.reserve(N);

    double vol = 1.0;
    for (double r : radii) vol *= r;
    if constexpr (Dim == 2) {
        vol *= M_PI;
    } else if constexpr (Dim == 3) {
        vol *= 4.0 / 3.0 * M_PI;
    }

    double per_particle = vol / N;
    double factor = (Dim == 2) ? 1.071 : std::pow(2.0, 1.0 / 6.0);
    double spacing = std::pow(per_particle, 1.0 / Dim) * factor;

    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-0.2 * spacing, 0.2 * spacing);

    int particles_added = 0;
    double buffer = spacing;

    if constexpr (Dim == 2) {
        double hex_offset = spacing * std::sqrt(3.0) / 2.0;
        int num_rows = static_cast<int>((2.0 * radii[1] + 2.0 * buffer) / hex_offset) + 2;

        for (int row = 0; row < num_rows && particles_added < N; ++row) {
            double y = center[1] - radii[1] - buffer + row * hex_offset;
            double x_offset = (row % 2 == 0) ? 0.0 : 0.5 * spacing;
            int num_cols = static_cast<int>((2.0 * radii[0] + 2.0 * buffer) / spacing) + 2;

            for (int col = 0; col < num_cols && particles_added < N; ++col) {
                double x = center[0] - radii[0] - buffer + col * spacing + x_offset;
                double dx = x - center[0];
                double dy = y - center[1];
                if (dx * dx / (radii[0] * radii[0]) + dy * dy / (radii[1] * radii[1]) <= 1.0) {
                    double jitter_x = dis(gen);
                    double jitter_y = dis(gen);
                    Particle<Dim> p;
                    p.pos[0] = x + jitter_x;
                    p.pos[1] = y + jitter_y;
                    particles.push_back(p);
                    particles_added++;
                }
            }
        }
    } else if constexpr (Dim == 3) {
        double hex_offset = spacing * std::sqrt(3.0) / 2.0;
        double layer_height = spacing * std::sqrt(6.0) / 3.0;
        int num_layers = static_cast<int>((2.0 * radii[2] + 2.0 * buffer) / layer_height) + 2;

        for (int layer = 0; layer < num_layers && particles_added < N; ++layer) {
            double z = center[2] - radii[2] - buffer + layer * layer_height;
            int layer_type = layer % 2;
            double x_layer_off = (layer_type == 1) ? spacing / 2.0 : 0.0;
            double y_layer_off = (layer_type == 1) ? spacing * std::sqrt(3.0) / 6.0 : 0.0;
            int num_rows = static_cast<int>((2.0 * radii[1] + 2.0 * buffer) / hex_offset) + 2;

            for (int row = 0; row < num_rows && particles_added < N; ++row) {
                double y = center[1] - radii[1] - buffer + row * hex_offset + y_layer_off;
                double x_row_off = (row % 2 == 0) ? 0.0 : spacing / 2.0;
                int num_cols = static_cast<int>((2.0 * radii[0] + 2.0 * buffer) / spacing) + 2;

                for (int col = 0; col < num_cols && particles_added < N; ++col) {
                    double x = center[0] - radii[0] - buffer + col * spacing + x_row_off + x_layer_off;
                    double dx = x - center[0];
                    double dy = y - center[1];
                    double dz = z - center[2];
                    if (dx * dx / (radii[0] * radii[0]) + dy * dy / (radii[1] * radii[1]) + dz * dz / (radii[2] * radii[2]) <= 1.0) {
                        double jitter_x = dis(gen);
                        double jitter_y = dis(gen);
                        double jitter_z = dis(gen);
                        Particle<Dim> p;
                        p.pos[0] = x + jitter_x;
                        p.pos[1] = y + jitter_y;
                        p.pos[2] = z + jitter_z;
                        particles.push_back(p);
                        particles_added++;
                    }
                }
            }
        }
    }

    if (particles_added < N) {
        std::cout << "Warning: Added only " << particles_added << " particles out of " << N << std::endl;
    }

    find_neighbors<Dim>(particles, h);
    compute_density<Dim>(particles, density0);

    double avg_rho = 0.0;
    for (const auto& p : particles) {
        avg_rho += p.density;
    }
    avg_rho /= particles.size();

    double scale = density0 / avg_rho;
    for (auto& p : particles) {
        p.mass *= scale;
    }
    compute_density<Dim>(particles, density0);
}

template void fill_ellipsoid<2>(Vector<2> center, std::array<double, 2> radii, std::vector<Particle<2>>& particles, int N, double density0, double h);
template void fill_ellipsoid<3>(Vector<3> center, std::array<double, 3> radii, std::vector<Particle<3>>& particles, int N, double density0, double h);