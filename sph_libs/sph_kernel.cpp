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

template<int Dim>
void initialize(double smoothing_length) {
    h = smoothing_length;
    h2 = h * h;
    h3 = h2 * h;
    inv_h = 1.0 / h;
    double h4 = h2 * h2;
    if (Dim == 2) {
        alpha_w = 7.0 / (4.0 * M_PI * h2);
        alpha_g = -35.0 / (4.0 * M_PI * h4);
    } else if (Dim == 3) {
        alpha_w = 21.0 / (16.0 * M_PI * h2);    // !Проверить!!!!!
        alpha_g = -105.0 / (16.0 * M_PI * h4);  // !Проверить!!!!!
    }

    h_CSS = 2 * smoothing_length;
    h2_CSS = h_CSS * h_CSS;
    h3_CSS = h2_CSS * h_CSS;
    inv_h_CSS = 1.0 / h_CSS;
    double h4_CSS = h2_CSS * h2_CSS;
    if (Dim == 2) {
        alpha_w_CSS = 7.0 / (4.0 * M_PI * h2_CSS);
        alpha_g_CSS = -35.0 / (4.0 * M_PI * h4_CSS);
    } else if (Dim == 3) {
        alpha_w_CSS = 21.0 / (16.0 * M_PI * h2_CSS);    // !Проверить!!!!!
        alpha_g_CSS = -105.0 / (16.0 * M_PI * h4_CSS);  // !Проверить!!!!!
    }
}

template<int Dim>
double wendlandC2(const Vector<Dim>& r) {
    double r_len = magnitude(r);
    if (r_len >= 2.0 * h) return 0.0;

    double q = r_len * inv_h;
    double tmp = 1.0 - 0.5 * q;
    double tmp2 = tmp * tmp;
    double tmp4 = tmp2 * tmp2;
    return alpha_w * tmp4 * (2.0 * q + 1.0);
}

template<int Dim>
Vector<Dim> grad_wendlandC2(const Vector<Dim>& r) {
    double r_len = magnitude(r);
    if (r_len == 0.0 || r_len >= 2.0 * h) return Vector<Dim>();

    double q = r_len * inv_h;
    double tmp = 1.0 - 0.5 * q;
    double tmp2 = tmp * tmp;
    double factor = tmp2 * tmp;
    double scale = alpha_g * factor;
    return r * scale;
}

template<int Dim>
double wendlandC2_CSS(const Vector<Dim>& r) {
    double r_len = magnitude(r);
    if (r_len == 0.0 || r_len >= 2.0 * h_CSS) return 0.0;

    double q = r_len * inv_h_CSS;
    double tmp = 1.0 - 0.5 * q;
    double tmp2 = tmp * tmp;
    double tmp4 = tmp2 * tmp2;
    return alpha_w_CSS * tmp4 * (2.0 * q + 1.0);
}

template<int Dim>
Vector<Dim> grad_wendlandC2_CSS(const Vector<Dim>& r) {
    double r_len = magnitude(r);
    if (r_len == 0.0 || r_len >= 2.0 * h_CSS) return Vector<Dim>();

    double q = r_len * inv_h_CSS;
    double tmp = 1.0 - 0.5 * q;
    double tmp2 = tmp * tmp;
    double factor = tmp2 * tmp;
    double scale = alpha_g_CSS * factor;
    return r * scale;
}

}  // namespace SPHKernel

template void SPHKernel::initialize<2>(double smoothing_length);
template void SPHKernel::initialize<3>(double smoothing_length);

template double SPHKernel::wendlandC2<2>(const Vector<2>& r);
template double SPHKernel::wendlandC2<3>(const Vector<3>& r);

template Vector<2> SPHKernel::grad_wendlandC2<2>(const Vector<2>& r);
template Vector<3> SPHKernel::grad_wendlandC2<3>(const Vector<3>& r);

template double SPHKernel::wendlandC2_CSS<2>(const Vector<2>& r);
template double SPHKernel::wendlandC2_CSS<3>(const Vector<3>& r);

template Vector<2> SPHKernel::grad_wendlandC2_CSS<2>(const Vector<2>& r);
template Vector<3> SPHKernel::grad_wendlandC2_CSS<3>(const Vector<3>& r);