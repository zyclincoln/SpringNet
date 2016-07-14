#include "point_potential_energy.cc"

using namespace Eigen;
using namespace std;
using namespace zyclincoln;

static const G = 0.001 * 9.8;

double PointPotentialEnergyCalculator::calculate_potential_energy(Point &point){
  return 0.5*point.mass_*pow(point.speed_, 2);
}

VectorXd PointPotentialEnergyCalculator::calculate_delta_potential_energy(Point &point){
  return Vector3d(0, -point.mass_*G, 0);
}

MatrixXd PointPotentialEnergyCalculator::calculate_delta_delta_potential_energy(Point &point){
  return MatrixXd::Zero(3, 3);
}