#include <sys/time.h>

#include <set>
#include <iostream>
#include <math.h>

#include "implicit_euler_integrator.h"

using namespace std;
using namespace Eigen;
using namespace zyclincoln;

void ImplicitEulerIntegrator::next_frame(AbstractSystem& system){
  timeval time0, time1, time2, time3, time4, time5, time6, time7;

  gettimeofday(&time0, 0);

  VectorXd f;
  system.delta_potential_energy_vector(f);

  gettimeofday(&time1, 0);
  MatrixXd stiff;
  system.delta_delta_potential_energy_matrix(stiff);

  gettimeofday(&time2, 0);
  if(solver_pre_computed_ == false || linear_ == false){
    MatrixXd A = system.mass_matrix() + 
      pow(system.time_step_ms(), 2)*stiff;
  /*
  
  vector<Triplet<double>> trips;
  trips.push_back(Triplet<double>(i, j, value));
  SparseMatrix<double> A;
  A.setFromTriplets(trips.begin(), trips.end());
  
  CSC format
  
  */
  
    // cout << (A-A.transpose()).norm() << endl;
    /* check symmetry */
    // cout << (A-A.transpose()).norm() << endl;
  
    solver.compute(A.sparseView());
    if(solver.info()!=Success){
      cout << "solver compute failed!" << endl;
      return;
    }
    solver_pre_computed_ = true;
    
    // cout <<"===A===\n" << A << endl;
  }

  // gettimeofday(&time0, 0);
  gettimeofday(&time3, 0);
  VectorXd b = -f*system.time_step_ms() - pow(system.time_step_ms(), 2)*stiff*system.velocity_vector();
  gettimeofday(&time4, 0);
  // cout <<"===f===\n" << system.delta_potential_energy_vector().transpose() << endl;
  VectorXd delta_velocity = solver.solve(b);
  // VectorXd delta_velocity = (system.mass_matrix() + pow(system.time_step_ms(), 2)*system.delta_delta_potential_energy_matrix()).inverse()*b;
  // gettimeofday(&time1,0);
  // double timeuse = 1000000*(time1.tv_sec - time0.tv_sec) + time1.tv_usec - time0.tv_usec;
  // cout << "solve time: " << timeuse/1000 << "ms" << endl;
  // cout << "f : \n" << f.transpose() << endl;
  // cout << "b : \n" << b.transpose() << endl;
  // cout << "delta velocity: \n" << delta_velocity.transpose() << endl;

  // system.total_energy();
  gettimeofday(&time5, 0);
  VectorXd new_velocity = system.velocity_vector() + delta_velocity;
  VectorXd new_position = system.position_vector() + new_velocity*system.time_step_ms();

  system.update_velocity_vector(new_velocity);
  system.update_position_vector(new_position);
  gettimeofday(&time6, 0);
  cout << "time to get f: " << (1000000*(time1.tv_sec - time0.tv_sec) + time1.tv_usec - time0.tv_usec)/1000 << " ms" << endl;
  cout << "time to get stiff: " << (1000000*(time2.tv_sec - time1.tv_sec) + time2.tv_usec - time1.tv_usec)/1000 << " ms" << endl;
  cout << "time to compute: " << (1000000*(time3.tv_sec - time2.tv_sec) + time3.tv_usec - time2.tv_usec)/1000 << " ms" << endl;
  cout << "time to get b: " << (1000000*(time4.tv_sec - time3.tv_sec) + time4.tv_usec - time3.tv_usec)/1000 << " ms" << endl;
  cout << "time to solve x: " << (1000000*(time5.tv_sec - time4.tv_sec) + time5.tv_usec - time4.tv_usec)/1000 << " ms" << endl;
  cout << "time to update: " << (1000000*(time6.tv_sec - time5.tv_sec) + time6.tv_usec - time5.tv_usec)/1000 << " ms" << endl;
  cout << "=====" << endl;
}