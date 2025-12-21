// neighbor_finder.cpp
#include "neighbor_finder.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "sph_kernel.hpp"

void find_neighbors(std::vector<particle>& particles, double h) {
    if (particles.empty()) return;

    // Clear existing neighbor lists
    for (auto& p : particles) {
        p.neighbors.clear();
    }

    double cell_size = 4.0 * h;
    double inv_cell_size = 1.0 / cell_size;
    double cutoff2 = 4.0 * h * h;

    // Find bounding box
    double min_x = particles[0].pos.x;
    double max_x = min_x;
    double min_y = particles[0].pos.y;
    double max_y = min_y;
    for (const auto& p : particles) {
        min_x = std::min(min_x, p.pos.x);
        max_x = std::max(max_x, p.pos.x);
        min_y = std::min(min_y, p.pos.y);
        max_y = std::max(max_y, p.pos.y);
    }

    // Add padding
    min_x -= cell_size;
    max_x += cell_size;
    min_y -= cell_size;
    max_y += cell_size;

    // Compute number of cells
    int num_cells_x = static_cast<int>(std::ceil((max_x - min_x) * inv_cell_size)) + 2;
    int num_cells_y = static_cast<int>(std::ceil((max_y - min_y) * inv_cell_size)) + 2;
    int total_cells = num_cells_x * num_cells_y;

    std::vector<std::vector<size_t>> cells(total_cells);

    // Assign particles to cells
    for (size_t i = 0; i < particles.size(); ++i) {
        double rel_x = particles[i].pos.x - min_x;
        double rel_y = particles[i].pos.y - min_y;
        int cx = static_cast<int>(rel_x * inv_cell_size);
        int cy = static_cast<int>(rel_y * inv_cell_size);
        if (cx >= 0 && cx < num_cells_x && cy >= 0 && cy < num_cells_y) {
            int cell_id = cy * num_cells_x + cx;
            cells[cell_id].push_back(i);
        }
    }

    // Build neighbor lists
    for (size_t i = 0; i < particles.size(); ++i) {
        double rel_x = particles[i].pos.x - min_x;
        double rel_y = particles[i].pos.y - min_y;
        int cx = static_cast<int>(rel_x * inv_cell_size);
        int cy = static_cast<int>(rel_y * inv_cell_size);

        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                int ncx = cx + dx;
                int ncy = cy + dy;
                if (ncx < 0 || ncx >= num_cells_x || ncy < 0 || ncy >= num_cells_y) continue;

                int ncell_id = ncy * num_cells_x + ncx;
                const auto& neigh_list = cells[ncell_id];

                for (size_t j : neigh_list) {
                    if (i == j) continue;

                    Vector2D rij = particles[i].pos - particles[j].pos;
                    double r2 = rij.x * rij.x + rij.y * rij.y;
                    if (r2 < 4 * cutoff2 && r2 > 0.0) {
                        particles[i].neighbors.push_back(j);
                    }
                }
            }
        }
    }
}