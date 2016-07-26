#include "point_potential_energy.h"

using namespace Eigen;
using namespace std;
using namespace zyclincoln;

//static const double G = 0.001 * 9.8;

static const double G = 0.05;

double PointPotentialEnergyCalculator::calculate_potential_energy(Point &point){
  return 0.5*point.mass_*pow(point.speed_.norm(), 2);
}

VectorXd PointPotentialEnergyCalculator::calculate_delta_potential_energy(Point &point){
  //return Vector3d(0, point.mass_*G, 0);
  return Vector3d::Zero();
}

MatrixXd PointPotentialEnergyCalculator::calculate_delta_delta_potential_energy(Point &point){
  return MatrixXd::Zero(3, 3);
}