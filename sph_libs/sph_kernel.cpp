// sph_kernel.cpp
#include "sph_kernel.hpp"

#include <cmath>

namespace SPHKernel {

double h = 0.1;
double h2 = 0.01;
double h3 = 0.001;
double inv_h = 10.0;
double alpha_w = 7.0 / (4.0 * M_PI * 0.01);
double alpha_g = -35.0 / (4.0 * M_PI * 0.01 * 0.01);

double h_CSS = 0.1;
double h2_CSS = 0.01;
double h3_CSS = 0.001;
double inv_h_CSS = 10.0;
double alpha_w_CSS = 7.0 / (4.0 * M_PI * 0.01);
double alpha_g_CSS = -35.0 / (4.0 * M_PI * 0.01 * 0.01);

void initialize(double smoothing_length) {
    h = smoothing_length;
    h2 = h * h;
    h3 = h2 * h;
    inv_h = 1.0 / h;
    double h4 = h2 * h2;
    alpha_w = 7.0 / (4.0 * M_PI * h2);
    alpha_g = -35.0 / (4.0 * M_PI * h4);

    h_CSS = 2 * smoothing_length;
    h2_CSS = h_CSS * h_CSS;
    h3_CSS = h2_CSS * h_CSS;
    inv_h_CSS = 1.0 / h_CSS;
    double h4_CSS = h2_CSS * h2_CSS;
    alpha_w_CSS = 7.0 / (4.0 * M_PI * h2_CSS);
    alpha_g_CSS = -35.0 / (4.0 * M_PI * h4_CSS);
}

double norm(const Vector2D& r) { return std::sqrt(r.x * r.x + r.y * r.y); }

double wendlandC2(const Vector2D& r) {
    double r2 = r.x * r.x + r.y * r.y;
    if (r2 >= 4.0 * h2) return 0.0;

    double q = std::sqrt(r2) * inv_h;
    double tmp = 1.0 - 0.5 * q;
    double tmp2 = tmp * tmp;
    double tmp4 = tmp2 * tmp2;

    return alpha_w * tmp4 * (2.0 * q + 1.0);
}

Vector2D grad_wendlandC2(const Vector2D& r) {
    double r2 = r.x * r.x + r.y * r.y;
    if (r2 >= 4.0 * h2 || r2 == 0.0) return {0.0, 0.0};

    double q = std::sqrt(r2) * inv_h;
    double tmp = 1.0 - 0.5 * q;
    double tmp2 = tmp * tmp;
    double factor = tmp2 * tmp;

    double scale = alpha_g * factor;
    return {scale * r.x, scale * r.y};
}

double wendlandC2_CSS(const Vector2D& r) {
    double r2 = r.x * r.x + r.y * r.y;
    if (r2 >= 4.0 * h2_CSS) return 0.0;

    double q = std::sqrt(r2) * inv_h_CSS;
    double tmp = 1.0 - 0.5 * q;
    double tmp2 = tmp * tmp;
    double tmp4 = tmp2 * tmp2;

    return alpha_w_CSS * tmp4 * (2.0 * q + 1.0);
}

Vector2D grad_wendlandC2_CSS(const Vector2D& r) {
    double r2 = r.x * r.x + r.y * r.y;
    if (r2 >= 4.0 * h2_CSS || r2 == 0.0) return {0.0, 0.0};

    double q = std::sqrt(r2) * inv_h_CSS;
    double tmp = 1.0 - 0.5 * q;
    double tmp2 = tmp * tmp;
    double factor = tmp2 * tmp;

    double scale = alpha_g_CSS * factor;
    return {scale * r.x, scale * r.y};
}

}  // namespace SPHKernel