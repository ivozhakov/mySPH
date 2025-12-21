// vtk_io.cpp
#include "vtk_io.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

void save_particles_to_vtk(const std::vector<particle>& particles, const std::string& filename, int timestep) {
    std::ofstream vtk_file(filename);
    if (!vtk_file.is_open()) {
        std::cerr << "Error: cannot open file " << filename << std::endl;
        return;
    }

    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << "Particle data at timestep " << timestep << "\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET UNSTRUCTURED_GRID\n";

    vtk_file << "POINTS " << particles.size() << " float\n";
    for (const auto& p : particles) {
        vtk_file << std::fixed << std::setprecision(10) << p.pos.x << " " << p.pos.y << " 0.0\n";
    }

    vtk_file << "CELLS " << particles.size() << " " << particles.size() * 2 << "\n";
    for (size_t i = 0; i < particles.size(); ++i) {
        vtk_file << "1 " << i << "\n";
    }

    vtk_file << "CELL_TYPES " << particles.size() << "\n";
    for (size_t i = 0; i < particles.size(); ++i) {
        vtk_file << "1\n";
    }

    vtk_file << "POINT_DATA " << particles.size() << "\n";

    vtk_file << "VECTORS velocity float\n";
    for (const auto& p : particles) {
        vtk_file << std::fixed << std::setprecision(10) << p.vel.x << " " << p.vel.y << " 0.0\n";
    }

    vtk_file << "SCALARS density float 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for (const auto& p : particles) {
        vtk_file << std::fixed << std::setprecision(10) << p.density << "\n";
    }

    vtk_file << "SCALARS pressure float 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for (const auto& p : particles) {
        vtk_file << std::fixed << std::setprecision(10) << p.pressure << "\n";
    }

    vtk_file << "SCALARS color float 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for (const auto& p : particles) {
        vtk_file << std::fixed << std::setprecision(10) << p.c << "\n";
    }

    vtk_file << "VECTORS grad_c float\n";
    for (const auto& p : particles) {
        vtk_file << std::fixed << std::setprecision(10) << p.grad_c.x << " " << p.grad_c.y << " 0.0\n";
    }

    vtk_file << "SCALARS kappa float 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for (const auto& p : particles) {
        vtk_file << std::fixed << std::setprecision(10) << p.kappa << "\n";
    }

    vtk_file.close();
}