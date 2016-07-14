#ifndef _POINT_POTENTIAL_ENERGY_H_
#define _POINT_POTENTIAL_ENERGY_H_

#include <Eigen/Dense>

#include "../model/point.h"

namespace zyclincoln{

  class PointPotentialEnergyCalculator{
  public:
    double calculate_potential_energy(Point &point);
    Eigen::VectorXd calculate_delta_potential_energy(Point &point);
    Eigen::MatrixXd calculate_delta_delta_potential_energy(Point &point);
  };
}

#endif