#include <Eigen/Sparse>
#include <set>
#include <math.h>

#include "implicit_euler_integrator.h"
#include "../utility/matrix_op.h"

using namespace std;
using namespace Eigen;
using namespace zyclincoln;

using index = unsigned int;

void ImplicitEulerIntegrator::next_frame(AbstractSystem& system){
  MatrixXd A = system.mass_matrix() + 
    pow(system.time_step_ms(), 2)*system.delta_delta_potential_energy_matrix();
  VectorXd b = system.delta_potential_energy_matrix() - 
    pow(system.time_step_ms(), 2)*system.delta_delta_potential_energy_matrix()*system.velocity_vector();

  VectorXd delta_velocity = A.colPivHouseholderQr().solve(b);

  VectorXd new_velocity = system.velocity_vector() + delta_velocity;
  VectorXd new_position = system.position_vector() + new_velocity*system.time_step_ms();

  system.update_velocity_vector(new_velocity);
  system.update_position_vector(new_positoin);
}