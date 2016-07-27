#include <sys/time.h>

#include <set>
#include <iostream>
#include <math.h>

#include "implicit_euler_integrator.h"

using namespace std;
using namespace Eigen;
using namespace zyclincoln;

void ImplicitEulerIntegrator::next_frame(AbstractSystem& system){
  timeval time0, time1;

  if(solver_pre_computed_ == false){
    MatrixXd A = system.mass_matrix() + 
      pow(system.time_step_ms(), 2)*system.delta_delta_potential_energy_matrix();
  /*
  
  vector<Triplet<double>> trips;
  trips.push_back(Triplet<double>(i, j, value));
  SparseMatrix<double> A;
  A.setFromTriplets(trips.begin(), trips.end());
  
  CSC format
  
  */
  
    // cout << (A-A.transpose()).norm() << endl;
    /* check symmetry */
    cout << (A-A.transpose()).norm() << endl;
  
    solver.compute(A.sparseView());
    if(solver.info()!=Success){
      cout << "solver compute failed!" << endl;
      return;
    }
    solver_pre_computed_ = true;
  }

  gettimeofday(&time0, 0);
  VectorXd b = -system.delta_potential_energy_vector()*system.time_step_ms() - 
  pow(system.time_step_ms(), 2)*system.delta_delta_potential_energy_matrix()*system.velocity_vector();

  VectorXd delta_velocity = solver.solve(b);

  gettimeofday(&time1,0);
  double timeuse = 1000000*(time1.tv_sec - time0.tv_sec) + time1.tv_usec - time0.tv_usec;
  cout << "solve time: " << timeuse/1000 << "ms" << endl;

  VectorXd new_velocity = system.velocity_vector() + delta_velocity;
  VectorXd new_position = system.position_vector() + new_velocity*system.time_step_ms();

  system.update_velocity_vector(new_velocity);
  system.update_position_vector(new_position);
}