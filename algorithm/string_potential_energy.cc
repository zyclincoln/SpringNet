#include "string_potential_energy.cc"

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;
using namespace zyclincoln;

double SpringPotentialEnergyCalculator::calculate_potential_energy(Spring &spring, Point &point0, Point &point1){
  return 0.5*spring.stiff_*pow((point0.position_ - point1.position_).norm(), 2);
}

VectorXd SpringPotentialEnergyCalculator::calculate_delta_potential_energy(Spring &spring, Point &point0, Point &point1){
  VectorXd delta = VectorXd::Zero(6, 1);
  delta.segment<3>(0) = -((point0.position_ - point1.position_).norm() - spring.length_)*spring.stiff_*(point0.position_ - point1.position_).normalized();
  delta.segment<3>(3) = -((point1.position_ - point0.position_).norm() - spring.length_)*spring.stiff_*(point1.position_ - point0.position_).normalized();

  return delta;
}

MatrixXd SpringPotentialEnergyCalculator::calculate_delta_delta_potential_energy(Spring &spring, Point &point0, Point &point1){
  MatrixXd delta = MatrixXd::Zero(6, 6);
  MatrixXd sub_delta = MatrixXd::Zero(3, 3);
  Vector3d p0 = point0.position_;
  Vector3d p1 = point1.position_;
  sub_delta = spring.stiff_*MatrixXd::Identity(3)*(1 - springs.length_/(p0 - p1).norm());
  sub_delta += spring.stiff_*spring.length_*(p1 - p0)*(p1 -p0).transpose()/pow((p1 - p0).norm(), 3);

  delta.block(0, 0, 3, 3) = sub_delta;
  delta.block(3, 3, 3, 3) = sub_delta;
  delta.block(0, 3, 3, 3) = -sub_delta;
  delta.block(3, 0, 3, 3) = -sub_delta;

  return delta;
}