#include "spring_potential_energy.h"

#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;
using namespace zyclincoln;

double SpringPotentialEnergyCalculator::calculate_potential_energy(Spring &spring, Point &point0, Point &point1){
  return 0.5*spring.stiff_*pow((point0.position_ - point1.position_).norm(), 2);
}

int SpringPotentialEnergyCalculator::calculate_delta_potential_energy(Spring &spring, Point &point0, Point &point1, VectorXd &delta){
  if(delta.rows() != 6 || delta.cols() != 1){
    delta.resize(6, 1);
  }
  delta = VectorXd::Zero(6, 1);
  delta.segment<3>(0) = ((point0.position_ - point1.position_).norm() - spring.length_)*spring.stiff_*(point0.position_ - point1.position_).normalized();
  delta.segment<3>(3) = ((point1.position_ - point0.position_).norm() - spring.length_)*spring.stiff_*(point1.position_ - point0.position_).normalized();

  return 0;
}

int SpringPotentialEnergyCalculator::calculate_delta_delta_potential_energy(Spring &spring, Point &point0, Point &point1, MatrixXd &delta){
  if(delta.rows() != 6 || delta.cols() != 6){
    delta.resize(6, 6);
  }
  delta = MatrixXd::Zero(6, 6);
  Matrix3d sub_delta = Matrix3d::Zero();
  Vector3d p0 = point0.position_;
  Vector3d p1 = point1.position_;
  sub_delta = spring.stiff_*MatrixXd::Identity(3, 3)*(1 - spring.length_/(p0 - p1).norm());
  sub_delta += spring.stiff_*spring.length_*(p1 - p0)*(p1 -p0).transpose()/pow((p1 - p0).norm(), 3);

  delta.block(0, 0, 3, 3) = sub_delta;
  delta.block(3, 3, 3, 3) = sub_delta;
  delta.block(0, 3, 3, 3) = -sub_delta;
  delta.block(3, 0, 3, 3) = -sub_delta;

  return 0;
}

double SimpleSpringPotentialEnergyCalculator::calculate_potential_energy(Spring &spring, Point &point0, Point &point1){
  return 0.5*spring.stiff_*pow((point0.position_ - point1.position_).norm(), 2);
}

int SimpleSpringPotentialEnergyCalculator::calculate_delta_potential_energy(Spring &spring, Point &point0, Point &point1, VectorXd &delta){
  if(delta.rows() != 3 || delta.cols() != 1){
    delta.resize(3, 1);
  }

  delta = spring.stiff_*(point0.position_ - point1.position_);

  return 0;
}

int SimpleSpringPotentialEnergyCalculator::calculate_delta_delta_potential_energy(Spring &spring, Point &point0, Point &point1, MatrixXd &delta){
  if(delta.rows() != 3 || delta.cols() != 3){
    delta.resize(3, 3);
  }

  delta = Matrix3d::Identity()*spring.stiff_;

  return 0;
}