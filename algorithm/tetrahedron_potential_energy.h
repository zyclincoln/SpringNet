#ifndef _TETRAHEDRON_POTENTIAL_ENERGY_H_
#define _TETRAHEDRON_POTENTIAL_ENERGY_H_

#include <Eigen/Dense>
#include "../model/tetrahedron.h"
#include "../model/point.h"

namespace zyclincoln{

  class TetrahedronPotentialEnergyCalculator{
  public:
    double calculate_potential_energy(Tetrahedron &tetrahedron);
    int calculate_delta_potential_energy(Tetrahedron &tetrahedron, Eigen::VectorXd &delta, bool corotate = false);
    int calculate_delta_delta_potential_energy(Tetrahedron &tetrahedron, Eigen::MatrixXd &delta, bool corotate = false);
    void PreCompute(Tetrahedron &tetrahedron, Point &point0, Point &point1, Point &point2, Point& point3, bool corotate = false);
  private:
    Eigen::MatrixXd delta_;
    double delta_0_;
  };
  
}

#endif