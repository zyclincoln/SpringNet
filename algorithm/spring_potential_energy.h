#ifndef _SPRING_POTENTIAL_ENERGY_H_
#define _SPRING_POTENTIAL_ENERGY_H_

#include <Eigen/Dense>

#include "../model/spring.h"
#include "../model/point.h"

namespace zyclincoln{

  class SpringPotentialEnergyCalculator{
  public:
    double calculate_potential_energy(Spring &spring, Point &point0, Point &point1);
    VectorXd calculate_delta_potential_energy(Spring &spring, Point &point0, Point &point1);
    MatrixXd calculator_delta_delta_potential_energy(Spring &spring, Point &point0, Point &point1);
  }
}

#endif