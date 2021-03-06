#ifndef _ABSTRACT_SYSTEM_H_
#define _ABSTRACT_SYSTEM_H_

#include <Eigen/Dense>

namespace zyclincoln{

  class AbstractSystem{
  public:
    virtual Eigen::MatrixXd& mass_matrix() = 0;
    virtual void delta_potential_energy_vector(Eigen::VectorXd &delta) = 0;
    virtual double total_energy() = 0;
    virtual void delta_delta_potential_energy_matrix(Eigen::MatrixXd &delta) = 0;
    virtual Eigen::VectorXd& velocity_vector() = 0;
    virtual Eigen::VectorXd& position_vector() = 0;
    virtual double time_step_ms() = 0;
    virtual void update_velocity_vector(const Eigen::VectorXd &velocity) = 0;
    virtual void update_position_vector(const Eigen::VectorXd &position) = 0;
  };

}

#endif