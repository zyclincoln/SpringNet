#ifndef _IMPLICIT_EULER_INTEGRATOR_H_
#define _IMPLICIT_EULER_INTEGRATOR_H_

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include "integrator.h"

namespace zyclincoln{

  class ImplicitEulerIntegrator : public Integrator{
  public:
    virtual void next_frame(AbstractSystem& system);
    void set_nonlinear();
  private:
    bool solver_pre_computed_ = false;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    bool linear_ = true;
  };

  inline void ImplicitEulerIntegrator::set_nonlinear(){
    linear_ = false;
  }

}

#endif