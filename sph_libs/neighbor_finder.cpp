// neighbor_finder.cpp
#include "neighbor_finder.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "sph_kernel.hpp"

template<int Dim>
int compute_cell_id(const std::array<int, Dim>& coords, const std::array<int, Dim>& num_cells) {
    int id = coords[0];
    int stride = num_cells[0];
    for (int d = 1; d < Dim; ++d) {
        id += coords[d] * stride;
        stride *= num_cells[d];
    }
    return id;
}

template<int Dim>
void find_neighbors(std::vector<Particle<Dim>>& particles, double h) {
    if (particles.empty()) return;

    // Clear existing neighbor lists
    for (auto& p : particles) {
        p.neighbors.clear();
    }

    double cutoff = 4.0 * h;
    double cutoff2 = cutoff * cutoff;
    double cell_size = cutoff;  // Cell size matches cutoff for efficiency
    double inv_cell_size = 1.0 / cell_size;

    // Find bounding box
    std::array<double, Dim> min_pos;
    std::array<double, Dim> max_pos;
    for (int d = 0; d < Dim; ++d) {
        min_pos[d] = particles[0].pos[d];
        max_pos[d] = min_pos[d];
    }
    for (const auto& p : particles) {
        for (int d = 0; d < Dim; ++d) {
            min_pos[d] = std::min(min_pos[d], p.pos[d]);
            max_pos[d] = std::max(max_pos[d], p.pos[d]);
        }
    }

    // Add padding
    for (int d = 0; d < Dim; ++d) {
        min_pos[d] -= cell_size;
        max_pos[d] += cell_size;
    }

    // Compute number of cells
    std::array<int, Dim> num_cells;
    int total_cells = 1;
    for (int d = 0; d < Dim; ++d) {
        num_cells[d] = static_cast<int>(std::ceil((max_pos[d] - min_pos[d]) * inv_cell_size)) + 2;
        total_cells *= num_cells[d];
    }

    std::vector<std::vector<size_t>> cells(total_cells);

    // Assign particles to cells
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& pos = particles[i].pos;
        std::array<int, Dim> cell_coords;
        bool valid = true;
        for (int d = 0; d < Dim; ++d) {
            double rel = pos[d] - min_pos[d];
            cell_coords[d] = static_cast<int>(rel * inv_cell_size);
            if (cell_coords[d] < 0 || cell_coords[d] >= num_cells[d]) {
                valid = false;
                break;
            }
        }
        if (valid) {
            int cell_id = compute_cell_id<Dim>(cell_coords, num_cells);
            cells[cell_id].push_back(i);
        }
    }

    // Compute number of offset combinations: 3^Dim
    int num_offsets = 1;
    for (int d = 0; d < Dim; ++d) {
        num_offsets *= 3;
    }

    // Build neighbor lists
    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& pi = particles[i];
        std::array<int, Dim> my_coords;
        for (int d = 0; d < Dim; ++d) {
            double rel = pi.pos[d] - min_pos[d];
            my_coords[d] = static_cast<int>(rel * inv_cell_size);
        }

        for (int off_id = 0; off_id < num_offsets; ++off_id) {
            int temp = off_id;
            std::array<int, Dim> offset;
            for (int d = 0; d < Dim; ++d) {
                offset[d] = (temp % 3) - 1;
                temp /= 3;
            }

            std::array<int, Dim> neigh_coords;
            bool valid = true;
            for (int d = 0; d < Dim; ++d) {
                neigh_coords[d] = my_coords[d] + offset[d];
                if (neigh_coords[d] < 0 || neigh_coords[d] >= num_cells[d]) {
                    valid = false;
                    break;
                }
            }
            if (!valid) continue;

            int neigh_cell_id = compute_cell_id<Dim>(neigh_coords, num_cells);
            const auto& neigh_list = cells[neigh_cell_id];

            for (size_t j : neigh_list) {
                if (i == j) continue;

                Vector<Dim> rij = pi.pos - particles[j].pos;
                double r2 = squared_norm(rij);
                
                if (r2 < cutoff2 && r2 > 1.e-15) {
                    particles[i].neighbors.push_back(j);
                }
            }
        }
    }
}

template void find_neighbors<2>(std::vector<Particle<2>>& particles, double h);
template void find_neighbors<3>(std::vector<Particle<3>>& particles, double h);