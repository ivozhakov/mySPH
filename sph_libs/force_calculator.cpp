// force_calculator.cpp
#include "force_calculator.hpp"

#include <cmath>
#include <vector>

#include "sph_kernel.hpp"

void compute_surface_tension_CSS(std::vector<particle>& particles, double sigma, double h) {
    if (sigma == 0.0 || particles.empty()) return;

    double w_self = SPHKernel::wendlandC2({0.0, 0.0});

    std::vector<Vector2D> normals(particles.size());

    // Compute grad_c with CSPM and KGC correction
    for (size_t i = 0; i < particles.size(); ++i) {
        const particle& pi = particles[i];
        double num_x = 0.0;
        double num_y = 0.0;
        double den = (pi.mass / pi.density) * w_self;

        double m00 = 0.0, m01 = 0.0, m10 = 0.0, m11 = 0.0;

        for (size_t j : pi.neighbors) {
            const particle& pj = particles[j];
            Vector2D rij = pi.pos - pj.pos;
            double r2 = rij.x * rij.x + rij.y * rij.y;
            if (r2 == 0.0) continue;

            double vol_j = pj.mass / pj.density;
            double w = SPHKernel::wendlandC2(rij);
            Vector2D gradW = SPHKernel::grad_wendlandC2(rij);

            num_x += vol_j * gradW.x;
            num_y += vol_j * gradW.y;
            den += vol_j * w;

            m00 += vol_j * rij.x * gradW.x;
            m01 += vol_j * rij.x * gradW.y;
            m10 += vol_j * rij.y * gradW.x;
            m11 += vol_j * rij.y * gradW.y;
        }

        double raw_gx = num_x;
        double raw_gy = num_y;

        double det = m00 * m11 - m01 * m10;
        if (std::abs(det) > 1e-10) {
            double inv_det = 1.0 / det;
            double l00 = m11 * inv_det;
            double l01 = -m01 * inv_det;
            double l10 = -m10 * inv_det;
            double l11 = m00 * inv_det;

            particles[i].grad_c.x = l00 * raw_gx + l01 * raw_gy;
            particles[i].grad_c.y = l10 * raw_gx + l11 * raw_gy;
        } else {
            if (den > 1e-6) {
                particles[i].grad_c.x = num_x / den;
                particles[i].grad_c.y = num_y / den;
            } else {
                particles[i].grad_c = Vector2D(0.0, 0.0);
            }
        }

        particles[i].cspm_den = den;
    }

    // Compute mag_grad_c and normals
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].mag_grad_c = SPHKernel::norm(particles[i].grad_c);
        if (particles[i].mag_grad_c < 1. / (4 * h)) {
            particles[i].grad_c = {0.0, 0.0};
            particles[i].mag_grad_c = 0.0;
        }

        if (particles[i].mag_grad_c > 1e-6) {
            normals[i] = particles[i].grad_c / particles[i].mag_grad_c;
        } else {
            normals[i] = Vector2D(0.0, 0.0);
        }
    }

    // Add symmetric surface tension forces using CSS approach
    for (size_t i = 0; i < particles.size(); ++i) {
        particle& pi = particles[i];
        for (size_t j : pi.neighbors) {
            const particle& pj = particles[j];
            Vector2D rij = pi.pos - pj.pos;
            double r2 = rij.x * rij.x + rij.y * rij.y;
            if (r2 == 0.0) continue;

            Vector2D gradW = SPHKernel::grad_wendlandC2_CSS(rij);

            double ndotg_i = normals[i].x * gradW.x + normals[i].y * gradW.y;
            Vector2D proj_i = gradW - normals[i] * ndotg_i;
            Vector2D dot_i = proj_i * (sigma * particles[i].mag_grad_c / (pi.density * pi.density));

            double ndotg_j = normals[j].x * gradW.x + normals[j].y * gradW.y;
            Vector2D proj_j = gradW - normals[j] * ndotg_j;
            Vector2D dot_j = proj_j * (sigma * particles[j].mag_grad_c / (pj.density * pj.density));

            Vector2D contrib = (dot_i + dot_j) * pj.mass;
            pi.force += contrib * pi.mass;
        }
    }
}

void density_normalization(std::vector<particle>& particles, double density0) {
    double w_self = SPHKernel::wendlandC2({0.0, 0.0});
    for (size_t i = 0; i < particles.size(); ++i) {
        const particle& pi = particles[i];
        double sum_vol_w = (pi.mass / pi.density) * w_self;
        for (size_t j : pi.neighbors) {
            Vector2D rij = pi.pos - particles[j].pos;
            double w = SPHKernel::wendlandC2(rij);
            sum_vol_w += (particles[j].mass / particles[j].density) * w;
        }
        if (sum_vol_w > 1e-6) {
            particles[i].density /= sum_vol_w;
        } else {
            particles[i].density = density0;
        }
    }
}

void compute_density(std::vector<particle>& particles, double density0) {
    double w_self = SPHKernel::wendlandC2({0.0, 0.0});
    for (size_t i = 0; i < particles.size(); ++i) {
        const particle& pi = particles[i];
        double rho = pi.mass * w_self;
        for (size_t j : pi.neighbors) {
            Vector2D rij = pi.pos - particles[j].pos;
            double w = SPHKernel::wendlandC2(rij);
            rho += particles[j].mass * w;
        }
        particles[i].density = std::max(rho, 1.e-6);
    }
    density_normalization(particles, density0);
}

void compute_pressure_forces(std::vector<particle>& particles, double cs, double density0, double h) {
    if (particles.empty()) return;

    double cutoff = 2.0 * h;
    double gamma = 7.0;
    double B = cs * cs * density0 / gamma;

    compute_density(particles, density0);

    for (auto& p : particles) {
        p.pressure = B * (std::pow(p.density / density0, gamma) - 1.0);
    }

    for (size_t i = 0; i < particles.size(); ++i) {
        const particle& pi = particles[i];
        for (size_t j : pi.neighbors) {
            if (j <= i) continue;

            const particle& pj = particles[j];
            Vector2D rij = pi.pos - pj.pos;
            double r2 = rij.x * rij.x + rij.y * rij.y;
            if (r2 == 0.0) continue;

            double r = std::sqrt(r2);
            if (r >= cutoff) continue;

            Vector2D gradW = SPHKernel::grad_wendlandC2(rij);

            double pressure_term = (pi.pressure / (pi.density * pi.density)) + (pj.pressure / (pj.density * pj.density));

            double mi_mj = pi.mass * pj.mass;
            Vector2D force_contrib = -(mi_mj * pressure_term) * gradW;

            particles[i].force += force_contrib;
            particles[j].force -= force_contrib;
        }
    }
}

void compute_viscous_forces(std::vector<particle>& particles, double viscosity, double h) {
    if (particles.empty()) return;

    double eta2 = 0.01 * h * h;
    double eight_visc = 8.0 * viscosity;

    for (size_t i = 0; i < particles.size(); ++i) {
        const particle& pi = particles[i];
        for (size_t j : pi.neighbors) {
            if (j <= i) continue;

            const particle& pj = particles[j];
            Vector2D rij = pi.pos - pj.pos;
            double r2 = rij.x * rij.x + rij.y * rij.y;
            if (r2 < 1.e-10) continue;

            // double r = std::sqrt(r2);
            Vector2D gradW = SPHKernel::grad_wendlandC2(rij);

            Vector2D vij = pi.vel - pj.vel;
            double vij_dot_rij = vij.x * rij.x + vij.y * rij.y;

            double visc_term = eight_visc * vij_dot_rij / (pi.density * pj.density * (r2 + eta2));

            double mi_mj = pi.mass * pj.mass;
            Vector2D force_contrib = mi_mj * visc_term * gradW;

            particles[i].force += force_contrib;
            particles[j].force -= force_contrib;
        }
    }
}