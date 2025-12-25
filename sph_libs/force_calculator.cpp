// force_calculator.cpp
#include "force_calculator.hpp"

#include <cmath>
#include <vector>
#include <optional>

#include "sph_kernel.hpp"

//внутренняя функция только для этого файла
template<int Dim>
std::optional<std::array<std::array<double, Dim>, Dim>>
invert_matrix(const std::array<std::array<double, Dim>, Dim>& m) {
    std::array<std::array<double, Dim>, Dim> inv{};
    double det = 0.0;
    
    if constexpr (Dim == 2) {
        // Детерминант для 2D
        det = m[0][0] * m[1][1] - m[0][1] * m[1][0];
        
        if (std::abs(det) <= 1e-10) {
            return std::nullopt; // Матрица вырождена
        }
        
        double inv_det = 1.0 / det;
        inv[0][0] =  m[1][1] * inv_det;
        inv[0][1] = -m[0][1] * inv_det;
        inv[1][0] = -m[1][0] * inv_det;
        inv[1][1] =  m[0][0] * inv_det;
        
    } else if constexpr (Dim == 3) {
        // Детерминант для 3D (исправленная формула)
        det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        
        if (std::abs(det) <= 1e-15) {
            return std::nullopt; // Матрица вырождена
        }
        
        double inv_det = 1.0 / det;
        
        // Вычисление алгебраических дополнений
        inv[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_det;
        inv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_det;
        inv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_det;
        
        inv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_det;
        inv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_det;
        inv[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_det;
        
        inv[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_det;
        inv[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_det;
        inv[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_det;
        
    } else {
        static_assert(Dim == 2 || Dim == 3, 
                     "Matrix inversion only implemented for 2D and 3D");
    }
    
    return inv;
}

template<int Dim>
void compute_surface_tension_CSS(std::vector<Particle<Dim>>& particles, double sigma, double h) {
    if (sigma == 0.0 || particles.empty() || h < 1.e-10) return;

    Vector<Dim> zero{};
    double w_self = SPHKernel::wendlandC2(zero);

    std::vector<Vector<Dim>> normals(particles.size());

    // Compute grad_c with CSPM and KGC correction
    for (size_t i = 0; i < particles.size(); ++i) {
        const Particle<Dim>& pi = particles[i];
        Vector<Dim> raw_grad{};
        //double num_x = 0.0;
        //double num_y = 0.0;
        double den = (pi.mass / pi.density) * w_self;

        std::array<std::array<double, Dim>, Dim> m{};
        //m[0][0] = 0.0;
        //m[1][0] = 0.0;
        //m[0][1] = 0.0;
        //m[1][1] = 0.0;

        for (size_t j : pi.neighbors) {
            const Particle<Dim>& pj = particles[j];
            Vector<Dim> rij = pi.pos - pj.pos;
            double r2 = squared_norm(rij);
            if (r2 < 1.e-15) continue;

            double vol_j = pj.mass / pj.density;
            double w = SPHKernel::wendlandC2(rij);
            Vector<Dim> gradW = SPHKernel::grad_wendlandC2(rij);

            raw_grad += vol_j * gradW;
            //num_x += vol_j * gradW[0];
            //num_y += vol_j * gradW[1];
            den += vol_j * w;

            for (int d = 0; d < Dim; ++d) {
                for (int e = 0; e < Dim; ++e) {
                    m[d][e] += vol_j * rij[d] * gradW[e];
                }
            }
            /*m[0][0] += vol_j * rij[0] * gradW[0];
            m[0][1] += vol_j * rij[0] * gradW[1];
            m[1][0] += vol_j * rij[1] * gradW[0];
            m[1][1] += vol_j * rij[1] * gradW[1];*/
        }

        //double raw_gx = num_x;
        //double raw_gy = num_y;

        std::optional<std::array<std::array<double, Dim>, Dim>> inv_m = invert_matrix<Dim>(m);

        if (inv_m.has_value()) {
            auto l = inv_m.value();
            particles[i].grad_c = Vector<Dim>{};
            // Матричное умножение l * raw_grad
            for (int d = 0; d < Dim; ++d) {
                for (int e = 0; e < Dim; ++e) {
                    particles[i].grad_c[d] += l[d][e] * raw_grad[e];
                }
            }
            //particles[i].grad_c[0] = l[0][0] * raw_gx + l[0][1] * raw_gy;
            //particles[i].grad_c[1] = l[1][0] * raw_gx + l[1][1] * raw_gy;
        } else {
            if (den > 1e-6) {
                particles[i].grad_c = raw_grad / den;
                //particles[i].grad_c[0] = num_x / den;
                //particles[i].grad_c[1] = num_y / den;
            } else {
                particles[i].grad_c = zero;
            }
        }

        /*double det = m00 * m11 - m01 * m10;
        if (std::abs(det) > 1e-10) {
            double inv_det = 1.0 / det;
            double l00 = m11 * inv_det;
            double l01 = -m01 * inv_det;
            double l10 = -m10 * inv_det;
            double l11 = m00 * inv_det;

            particles[i].grad_c[0] = l00 * raw_gx + l01 * raw_gy;
            particles[i].grad_c[1] = l10 * raw_gx + l11 * raw_gy;
        } else {
            if (den > 1e-6) {
                particles[i].grad_c[0] = num_x / den;
                particles[i].grad_c[1] = num_y / den;
            } else {
                particles[i].grad_c = zero;
            }
        }*/

        particles[i].cspm_den = den;
    }

    // Compute mag_grad_c and normals
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].mag_grad_c = magnitude(particles[i].grad_c);
        if (particles[i].mag_grad_c < 0.5*SPHKernel::inv_h_CSS) {
            particles[i].grad_c = zero;
            particles[i].mag_grad_c = 0.0;
        }

        if (particles[i].mag_grad_c > 1e-6) {
            normals[i] = particles[i].grad_c / particles[i].mag_grad_c;
        } else {
            normals[i] = zero;
        }
    }

    // Add symmetric surface tension forces using CSS approach
    for (size_t i = 0; i < particles.size(); ++i) {
        Particle<Dim>& pi = particles[i];
        for (size_t j : pi.neighbors) {
            const Particle<Dim>& pj = particles[j];
            Vector<Dim> rij = pi.pos - pj.pos;
            double r2 = squared_norm(rij);
            if (r2 == 0.0) continue;

            Vector<Dim> gradW = SPHKernel::grad_wendlandC2_CSS(rij);

            double ndotg_i = dot(normals[i], gradW);
            //double ndotg_i = normals[i][0] * gradW[0] + normals[i][1] * gradW[1];
            Vector<Dim> proj_i = gradW - normals[i] * ndotg_i;
            Vector<Dim> dot_i = proj_i * (sigma * particles[i].mag_grad_c / (pi.density * pi.density));

            double ndotg_j = dot(normals[j], gradW);
            //double ndotg_j = normals[j][0] * gradW[0] + normals[j][1] * gradW[1];
            Vector<Dim> proj_j = gradW - normals[j] * ndotg_j;
            Vector<Dim> dot_j = proj_j * (sigma * particles[j].mag_grad_c / (pj.density * pj.density));

            Vector<Dim> contrib = (dot_i + dot_j) * pj.mass;
            pi.force += contrib * pi.mass;
        }
    }
}

template<int Dim>
void density_normalization(std::vector<Particle<Dim>>& particles, double density0) {
    Vector<Dim> zero{};
    double w_self = SPHKernel::wendlandC2<Dim>(zero);
    for (size_t i = 0; i < particles.size(); ++i) {
        Particle<Dim>& pi = particles[i];
        double sum_vol_w = (pi.mass / pi.density) * w_self;
        for (size_t j : pi.neighbors) {
            Particle<Dim>& pj = particles[j];
            Vector<Dim> rij = pi.pos - pj.pos;
            double w = SPHKernel::wendlandC2<Dim>(rij);
            sum_vol_w += (pj.mass / pj.density) * w;
        }
        if (sum_vol_w > 1e-6) {
            pi.density /= sum_vol_w;
        } else {
            pi.density = density0;
        }
    }
}

template<int Dim>
void compute_density(std::vector<Particle<Dim>>& particles, double density0) {
    Vector<Dim> zero{};
    double w_self = SPHKernel::wendlandC2<Dim>(zero);
    for (size_t i = 0; i < particles.size(); ++i) {
        Particle<Dim>& pi = particles[i];
        double rho = pi.mass * w_self;
        for (size_t j : pi.neighbors) {
            Particle<Dim>& pj = particles[j];
            Vector<Dim> rij = pi.pos - pj.pos;
            double w = SPHKernel::wendlandC2<Dim>(rij);
            rho += pj.mass * w;
        }
        pi.density = std::max(rho, 1.e-6);
    }
    density_normalization(particles, density0);
}

template<int Dim>
void compute_pressure_forces(std::vector<Particle<Dim>>& particles, double cs, double density0, double h) {
    if (particles.empty()) return;

    double cutoff = 2.0 * h;
    double gamma = 7.0;
    double B = cs * cs * density0 / gamma;

    compute_density<Dim>(particles, density0);

    for (auto& p : particles) {
        p.pressure = B * (std::pow(p.density / density0, gamma) - 1.0);
    }

    for (size_t i = 0; i < particles.size(); ++i) {
        const Particle<Dim>& pi = particles[i];
        for (size_t j : pi.neighbors) {
            if (j <= i) continue;

            const Particle<Dim>& pj = particles[j];
            Vector<Dim> rij = pi.pos - pj.pos;
            double r = magnitude(rij);
            if (r >= cutoff || r == 0.0) continue;

            Vector<Dim> gradW = SPHKernel::grad_wendlandC2<Dim>(rij);

            double pressure_term = (pi.pressure / (pi.density * pi.density)) + (pj.pressure / (pj.density * pj.density));

            double mi_mj = pi.mass * pj.mass;
            Vector<Dim> force_contrib = -(mi_mj * pressure_term) * gradW;

            particles[i].force += force_contrib;
            particles[j].force -= force_contrib;
        }
    }
}

template<int Dim>
void compute_viscous_forces(std::vector<Particle<Dim>>& particles, double viscosity, double h) {
    if (particles.empty() || h == 0.0) return;

    double eta2 = 0.01 * SPHKernel::h2;
    double visc_coeff = 2.0 * (static_cast<double>(Dim) + 2.0) * viscosity;  // 8 for 2D, 10 for 3D

    for (size_t i = 0; i < particles.size(); ++i) {
        const Particle<Dim>& pi = particles[i];
        for (size_t j : pi.neighbors) {
            if (j <= i) continue;

            const Particle<Dim>& pj = particles[j];
            Vector<Dim> rij = pi.pos - pj.pos;
            double r2 = squared_norm(rij);
            if (r2 < 1.e-10) continue;

            // double r = std::sqrt(r2);
            Vector<Dim> gradW = SPHKernel::grad_wendlandC2<Dim>(rij);

            Vector<Dim> vij = pi.vel - pj.vel;
            double vij_dot_rij = dot(vij, rij);

            double visc_term = visc_coeff * vij_dot_rij / (pi.density * pj.density * (r2 + eta2));

            double mi_mj = pi.mass * pj.mass;
            Vector<Dim> force_contrib = mi_mj * visc_term * gradW;

            particles[i].force += force_contrib;
            particles[j].force -= force_contrib;
        }
    }
}

template void compute_surface_tension_CSS<2>(std::vector<Particle<2>>& particles, double sigma, double h);
template void compute_surface_tension_CSS<3>(std::vector<Particle<3>>& particles, double sigma, double h);

template void compute_density<2>(std::vector<Particle<2>>& particles, double density0);
template void compute_density<3>(std::vector<Particle<3>>& particles, double density0);

template void density_normalization<2>(std::vector<Particle<2>>& particles, double density0);
template void density_normalization<3>(std::vector<Particle<3>>& particles, double density0);

template void compute_pressure_forces<2>(std::vector<Particle<2>>& particles, double cs, double density0, double h);
template void compute_pressure_forces<3>(std::vector<Particle<3>>& particles, double cs, double density0, double h);

template void compute_viscous_forces<2>(std::vector<Particle<2>>& particles, double viscosity, double h);
template void compute_viscous_forces<3>(std::vector<Particle<3>>& particles, double viscosity, double h);