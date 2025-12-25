// sph_solver.hpp
#ifndef SPH_SOLVER_HPP
#define SPH_SOLVER_HPP

#include <vector>

#include "particle.hpp"
#include "vector.hpp"

template<int Dim>
class SPHSolver {
   private:
    double h;
    double dt;
    double density0;
    double cs;
    double viscosity;
    double sigma;
    Vector<Dim> gravity;
    double box_length;

   public:
    SPHSolver(double smoothing_length = 0.1, double time_step = 0.00001, double ref_density = 1000.0, double cs = 10.0, double fluid_viscosity = 0.0, double surface_tension = 0.07,
              Vector<Dim> acceleration = Vector<Dim>(), double box_size = 0.005);

    void compute_forces(std::vector<Particle<Dim>>& particles);
    void simulate_step(std::vector<Particle<Dim>>& particles);

    double get_smoothing_length() const;
    double get_time_step() const;

   private:
    void apply_boundaries(std::vector<Particle<Dim>>& particles) const;
};

#endif  // SPH_SOLVER_HPP