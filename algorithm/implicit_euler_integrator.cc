 /*
 
 vector<Triplet<double>> trips;
 trips.push_back(Triplet<double>(i, j, value));
 SparseMatrix<double> A;
 A.setFromTriplets(trips.begin(), trips.end());
 
 CSC format
 
 */

#include <sys/time.h>

#include <set>
#include <iostream>
#include <math.h>
#include <Eigen/Eigenvalues>

#include "implicit_euler_integrator.h"
#include "../utility/eigen_util.h"

using namespace std;
using namespace Eigen;
using namespace zyclincoln;

void ImplicitEulerIntegrator::next_frame(AbstractSystem& system){
  if(reduce_ == false){
    next_frame_without_reduce(system);
  }
  else{
    next_frame_with_reduce(system);
  }
}

void ImplicitEulerIntegrator::next_frame_without_reduce(AbstractSystem& system){
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

    solver.compute(A.sparseView());
    if(solver.info()!=Success){
      cout << "solver compute failed!" << endl;
      return;
    }
    solver_pre_computed_ = true;
  }

  gettimeofday(&time3, 0);
  VectorXd b = -f*system.time_step_ms() - pow(system.time_step_ms(), 2)*stiff*system.velocity_vector();
  gettimeofday(&time4, 0);
  VectorXd delta_velocity = solver.solve(b);

  gettimeofday(&time5, 0);
  VectorXd new_velocity = system.velocity_vector() + delta_velocity;
  VectorXd new_position = system.position_vector() + new_velocity*system.time_step_ms();

  system.update_velocity_vector(new_velocity);
  system.update_position_vector(new_position);
  gettimeofday(&time6, 0);
  // cout << "time to get f: " << (1000000*(time1.tv_sec - time0.tv_sec) + time1.tv_usec - time0.tv_usec)/1000 << " ms" << endl;
  // cout << "time to get stiff: " << (1000000*(time2.tv_sec - time1.tv_sec) + time2.tv_usec - time1.tv_usec)/1000 << " ms" << endl;
  // cout << "time to compute: " << (1000000*(time3.tv_sec - time2.tv_sec) + time3.tv_usec - time2.tv_usec)/1000 << " ms" << endl;
  // cout << "time to get b: " << (1000000*(time4.tv_sec - time3.tv_sec) + time4.tv_usec - time3.tv_usec)/1000 << " ms" << endl;
  // cout << "time to solve x: " << (1000000*(time5.tv_sec - time4.tv_sec) + time5.tv_usec - time4.tv_usec)/1000 << " ms" << endl;
  // cout << "time to update: " << (1000000*(time6.tv_sec - time5.tv_sec) + time6.tv_usec - time5.tv_usec)/1000 << " ms" << endl;
  // cout << "=====" << endl;
}

void ImplicitEulerIntegrator::next_frame_with_reduce(AbstractSystem& system){
  VectorXd f;
  system.delta_potential_energy_vector(f);

  MatrixXd stiff;
  system.delta_delta_potential_energy_matrix(stiff);

  if(solver_pre_computed_ == false){
    EigenSolver<MatrixXd> eigen_solver;
    cout << "===begin solve eigen===\n" << endl;
    eigen_solver.compute((system.mass_matrix().inverse()*stiff), true);
    cout << "===finish===\n" << endl;
    VectorXcd eigen_complex_value = eigen_solver.eigenvalues();
    MatrixXcd eigen_complex_vector = eigen_solver.eigenvectors();
    VectorXd eigen_value;
    MatrixXd eigen_vector;

    eigen_value.resize(eigen_complex_value.rows(), eigen_complex_value.cols());
    eigen_vector.resize(eigen_complex_vector.rows(), eigen_complex_vector.cols());

    for(unsigned int i = 0; i < eigen_value.cols(); i++){
      eigen_value(i, 0) = eigen_complex_value(i, 0).real();
      if(fabs(eigen_complex_value(i, 0).imag()) > 10e-9){
        cerr << "[ERROR] eigen value is not real." << endl;
      }
      for(unsigned int j = 0 ; j < eigen_complex_vector.rows(); j++){
        eigen_vector(j, i) = eigen_complex_vector(j, i).real();
        if(fabs(eigen_complex_vector(j, i).imag()) > 10e-9){
          cerr << "[ERROR] eigen vector is not real." << endl;
        }
      }
    }

    cout << "value:\n " << eigen_value << endl;
    cout << "vector:\n " << eigen_vector << endl;
    
    get_top_k_eigen(eigen_value, eigen_vector, reduce_eigen_value_, reduce_eigen_vector_, 100);

    MatrixXd A = reduce_eigen_vector_.transpose()*system.mass_matrix()*reduce_eigen_vector_ + reduce_eigen_vector_.transpose()*stiff*reduce_eigen_vector_*pow(system.time_step_ms(), 2);
    solver.compute(A.sparseView());
    if(solver.info() != Success){
      cout << "solver compute failed!" << endl;
      return;
    }
    solver_pre_computed_ = true;
  }

  // VectorXd b = -reduce_eigen_vector_.transpose()*f*system.time_step_ms() - pow(system.time_step_ms(), 2)*(reduce_eigen_vector_.transpose()*stiff*reduce_eigen_vector_)*(reduce_eigen_vector_*system.velocity_vector());
  // VectorXd delta_velocity = reduce_eigen_vector_.transpose()*solver.solve(b);

  // VectorXd new_velocity = system.velocity_vector() + delta_velocity;
  // VectorXd new_position = system.position_vector() + new_velocity*system.time_step_ms();

  // system.update_velocity_vector(new_velocity);
  // system.update_position_vector(new_position);
}