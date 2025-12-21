// sph_solver.hpp
#ifndef SPH_SOLVER_HPP
#define SPH_SOLVER_HPP

#include <vector>

#include "particle.hpp"
#include "vector2d.hpp"

class SPHSolver {
   private:
    double h;
    double dt;
    double density0;
    double stiffness;
    double viscosity;
    double sigma;
    Vector2D gravity;
    double box_length;

   public:
    SPHSolver(double smoothing_length = 0.1, double time_step = 0.00001, double ref_density = 1000.0, double pressure_stiffness = 1000.0, double fluid_viscosity = 0.0, double surface_tension = 0.07,
              Vector2D acceleration = {0.0, 0.0}, double box_size = 0.005);

    void compute_forces(std::vector<particle>& particles);
    void simulate_step(std::vector<particle>& particles);

    double get_smoothing_length() const;
    double get_time_step() const;

   private:
    void apply_boundaries(std::vector<particle>& particles) const;
};

#endif  // SPH_SOLVER_HPP