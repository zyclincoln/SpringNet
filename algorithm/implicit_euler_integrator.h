#ifndef _IMPLICIT_EULER_INTEGRATOR_H_
#define _IMPLICIT_EULER_INTEGRATOR_H_

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include "integrator.h"

namespace zyclincoln{

  class ImplicitEulerIntegrator : public Integrator{
  public:
    virtual void next_frame(AbstractSystem& system);
    void next_frame_without_reduce(AbstractSystem& system);
    void next_frame_with_reduce(AbstractSystem& system);
    void set_nonlinear();
    void set_reduce();
  private:
    bool solver_pre_computed_ = false;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    bool linear_ = true;
    bool reduce_ = false;
    Eigen::MatrixXd reduce_eigen_vector_;
    Eigen::VectorXd reduce_eigen_value_;
  };

  inline void ImplicitEulerIntegrator::set_nonlinear(){
    linear_ = false;
  }

  inline void ImplicitEulerIntegrator::set_reduce(){
    reduce_ = true;
  }

}

#endif