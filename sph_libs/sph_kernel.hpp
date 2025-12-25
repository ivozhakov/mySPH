// sph_kernel.hpp
#ifndef SPH_KERNEL_HPP
#define SPH_KERNEL_HPP

#include "vector.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace SPHKernel {

extern double h;
extern double h2;
extern double h3;
extern double inv_h;
extern double alpha_w;
extern double alpha_g;

extern double h_CSS;
extern double h2_CSS;
extern double h3_CSS;
extern double inv_h_CSS;
extern double alpha_w_CSS;
extern double alpha_g_CSS;

template<int Dim>
void initialize(double smoothing_length);

template<int Dim>
double wendlandC2(const Vector<Dim>& r);

template<int Dim>
Vector<Dim> grad_wendlandC2(const Vector<Dim>& r);

template<int Dim>
double wendlandC2_CSS(const Vector<Dim>& r);

template<int Dim>
Vector<Dim> grad_wendlandC2_CSS(const Vector<Dim>& r);

}  // namespace SPHKernel

#endif  // SPH_KERNEL_HPP