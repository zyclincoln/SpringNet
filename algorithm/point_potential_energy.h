#ifndef _POINT_POTENTIAL_ENERGY_H_
#define _POINT_POTENTIAL_ENERGY_H_

#include <Eigen/Dense>

namespace zyclincoln{

  class PointPotentialEnergyCalculator{
  public:
    double claculate_potential_energy(Point &point);
    VectorXd calculate_delta_potential_energy(Point &point);
    MatrixXd calculator_delta_delta_potential_energy(Point &point);
  }
}

#endif