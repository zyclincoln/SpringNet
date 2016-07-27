/*
X(3,4)
dx(col(x,2)-col(x,1), col(x,3)-col(x,1), col(x,4)-col(x,1));

d=dx^(-1)
v=(-d11-d21-d31, -d12-d22-d32, -d13-d23-d33)

d2P = {
sub(1,1) = (2*lambda+miu)*v^t*v + miu*v*v^t*I(3);
sub(1,2) = miu*(d11, d12, d13)^t*v + 2*lambda*v^t*(d11, d12, d13) + miu*(d11, d12, d13)*v^t*I(3);
sub(1,3) = miu*(d21, d22, d23)^t*v + 2*lambda*v^t*(d21, d22, d23) + miu*(d21, d22, d23)*v^t*I(3);
sub(1,4) = miu*(d31, d32, d33)^t*v + 2*lambda*v^t*(d31, d32, d33) + miu*(d31, d32, d33)*v^t*I(3);
sub(2,1) = sub(1,2)^t;
sub(2,2) = 2*lambda*(d11, d12, d13)^t*(d11, d12, d13) + miu*(d11, d12, d13)^t*(d11, d12, d13) + miu*(d11, d12, d13)*(d11, d12, d13)^t*I(3);
sub(2,3) = 2*lambda*(d11, d12, d13)^t*(d21, d22, d23) + miu*(d21, d22, d23)^t*(d11, d12, d13) + miu*(d11, d12, d13)*(d21, d22, d23)^t*I(3);
sub(2,4) = 2*lambda*(d11, d12, d13)^t*(d31, d32, d33) + miu*(d31, d32, d33)^t*(d11, d12, d13) + miu*(d11, d12, d13)*(d31, d32, d33)^t*I(3);
sub(3,2) = sub(2,3)^t;
}
*/

/*

dP = {
dP(1,1) = 2*lambda*v1*(-3 + delta0) + miu*(2*v1*delta(1,1) + v2*delta(1,2) + v3*delta(1,3));
dP(2,1) = 2*lambda*v2*(-3 + delta0) + miu*(v1*delta(2,1) + 2*v2*delta(2,2) + v3*delta(2,3));
dP(3,1) = 2*lambda*v3*(-3 + delta0) + miu*(v1*delta(3,1) + v2*delta(3,2) + 2*v3*delta(3,3));

dP(4,1) = 2*lambda*d11*(-3 + delta0) + miu*(2*d11*delta(1,1) + d12*delta(1,2) + d13*delta(1,3));
dP(5,1) = 2*lambda*d12*(-3 + delta0) + miu*(d11*delta(2,1) + 2*d12*delta(2,2) + d13*delta(2,3));
...
dP(9,1) = 2*lambda*d23*(-3 + delta0) + miu*(d21*delta(3,1) + d22*delta(3,2) + 2*d23*delta(3,3));
}

where

delta0 = (d11 d21 d31 d12 d22 d32 d13 d23 d33)*(dx11 dx12 dx13 dx21 dx22 dx23 dx31 dx32 dx33)^t;

delta(1,1) = d31*dx13 + d21*dx12 + d11*dx11 - 1;
delta(2,2) = d32*dx23 + d22*dx22 + d12*dx21 - 1;
delta(3,3) = d33*dx33 + d23*dx32 + d13*dx31 - 1;
delta(1,2) = delta(2,1) = d32*dx13 + d31*dx23 + d21*dx22 + d22*dx12 + d11*dx21 + d12*dx11;
delta(1,3) = delta(3,1) = d31*dx33 + d21*dx32 + d33*dx13 + d11*dx31 + d23*dx12 + d13*dx11;
delta(2,3) = delta(3,2) = d32*dx33 + d22*dx32 + d33*dx23 + d23*dx22 + d12*dx31 + d13*dx21;


P = lambda*(-3+delta0)^2 + miu*(delta(1,1)^2 + delta(1,2)^2/2 + delta(2,2)^2 + delta(1,3)^2/2 + delta(3,3)^2 + delta(2,3)^2/2);

*/

#include "tetrahedron_potential_energy.h"

#include <iostream>
#include <math.h>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
using namespace zyclincoln;

void TetrahedronPotentialEnergyCalculator::PreCompute(Tetrahedron &tetrahedron, Point &point0, Point &point1, Point &point2, Point& point3){
  vector<Point> points_set;
  points_set.push_back(point0);
  points_set.push_back(point1);
  points_set.push_back(point2);
  points_set.push_back(point3);
  tetrahedron.ComputeDx(points_set);
  // cout << "inverse check : " << endl << tetrahedron.dx_*tetrahedron.p_inverse_ << endl;
  delta_ = Matrix3d::Zero();

  delta_(0, 0) = tetrahedron.p_inverse_(2, 0)*tetrahedron.dx_(0, 2) + tetrahedron.p_inverse_(1, 0)*tetrahedron.dx_(0, 1)
    + tetrahedron.p_inverse_(0, 0)*tetrahedron.dx_(0, 0) - 1;
  delta_(1, 1) = tetrahedron.p_inverse_(2, 1)*tetrahedron.dx_(1, 2) + tetrahedron.p_inverse_(1, 1)*tetrahedron.dx_(1, 1)
    + tetrahedron.p_inverse_(0, 1)*tetrahedron.dx_(1, 0) - 1;
  delta_(2, 2) = tetrahedron.p_inverse_(2, 2)*tetrahedron.dx_(2, 2) + tetrahedron.p_inverse_(1, 2)*tetrahedron.dx_(2, 1)
    + tetrahedron.p_inverse_(0, 2)*tetrahedron.dx_(2, 0) - 1;
  delta_(0, 1) = tetrahedron.p_inverse_(2, 1)*tetrahedron.dx_(0, 2) + tetrahedron.p_inverse_(2, 0)*tetrahedron.dx_(1, 2)
    + tetrahedron.p_inverse_(1, 0)*tetrahedron.dx_(1, 1) + tetrahedron.p_inverse_(1, 1)*tetrahedron.dx_(0, 1)
    + tetrahedron.p_inverse_(0, 0)*tetrahedron.dx_(1, 0) + tetrahedron.p_inverse_(0, 1)*tetrahedron.dx_(0, 0);
  delta_(0, 2) = tetrahedron.p_inverse_(2, 0)*tetrahedron.dx_(2, 2) + tetrahedron.p_inverse_(2, 2)*tetrahedron.dx_(0, 2)
    + tetrahedron.p_inverse_(1, 0)*tetrahedron.dx_(2, 1) + tetrahedron.p_inverse_(1, 2)*tetrahedron.dx_(0, 1)
    + tetrahedron.p_inverse_(0, 0)*tetrahedron.dx_(2, 0) + tetrahedron.p_inverse_(0, 2)*tetrahedron.dx_(0, 0);
  delta_(1, 2) = tetrahedron.p_inverse_(2, 1)*tetrahedron.dx_(2, 2) + tetrahedron.p_inverse_(2, 2)*tetrahedron.dx_(1, 2)
    + tetrahedron.p_inverse_(1, 1)*tetrahedron.dx_(2, 1) + tetrahedron.p_inverse_(1, 2)*tetrahedron.dx_(1, 1)
    + tetrahedron.p_inverse_(0, 1)*tetrahedron.dx_(2, 0) + tetrahedron.p_inverse_(0, 2)*tetrahedron.dx_(1, 0);
  delta_(2, 1) = delta_(1, 2);
  delta_(2, 0) = delta_(0, 2);
  delta_(1, 0) = delta_(0, 1);

  // delta_0_ =  tetrahedron.p_inverse_(0, 0)*tetrahedron.dx_(0, 0) + tetrahedron.p_inverse_(1, 0)*tetrahedron.dx_(0, 1)
  //   + tetrahedron.p_inverse_(2, 0)*tetrahedron.dx_(0, 2) + tetrahedron.p_inverse_(0, 1) + tetrahedron.dx_(1, 0)
  //   + tetrahedron.p_inverse_(1, 1)*tetrahedron.dx_(1, 1) + tetrahedron.p_inverse_(2, 1) + tetrahedron.dx_(1, 2)
  //   + tetrahedron.p_inverse_(0, 2)*tetrahedron.dx_(2, 0) + tetrahedron.p_inverse_(1, 2) + tetrahedron.dx_(2, 1)
  //   + tetrahedron.p_inverse_(2, 2)*tetrahedron.dx_(2, 2);

  // delta_0_ = (tetrahedron.p_inverse_.row(0)*tetrahedron.dx_.col(0) + tetrahedron.p_inverse_.row(1)*tetrahedron.dx_.col(1) + tetrahedron.p_inverse_.row(2)*tetrahedron.dx_.col(2))(0,0);
  delta_0_ = (tetrahedron.p_inverse_*tetrahedron.dx_).trace();
  // cout << "delta_0_: " << delta_0_ << endl;
}

double TetrahedronPotentialEnergyCalculator::calculate_potential_energy(Tetrahedron &tetrahedron){
  double energy = 0;
  energy = tetrahedron.lambda_*pow((-3 + delta_0_), 2);
  energy += tetrahedron.miu_*(pow(delta_(0, 0), 2) + pow(delta_(0, 1), 2)*0.5 + pow(delta_(1, 1), 2)*0.5
    + pow(delta_(0, 2), 2)*0.5 + pow(delta_(2, 2), 2) + pow(delta_(1, 2), 2)*0.5);
  // cout << "energy : " << energy << endl;
  // cout << "delta_0_: " << delta_0_ << endl;
  // cout << "delta 0,1: " << delta_(0, 1) << endl;
  // cout << "delta 0,2: " << delta_(0, 2) << endl;
  // cout << "delta 1,2: " << delta_(1, 2) << endl;
  // cout << "delta 0,0: " << delta_(0, 0) << endl;
  // cout << "delta 1,1: " << delta_(1, 1) << endl;
  // cout << "delta 2,2: " << delta_(2, 2) << endl;
  // cout << "miu : " << tetrahedron.miu_ << endl;
  // cout << "lambda : " << tetrahedron.lambda_ << endl;
  return energy;
}

VectorXd TetrahedronPotentialEnergyCalculator::calculate_delta_potential_energy(Tetrahedron &tetrahedron){
  VectorXd elastic_force;
  elastic_force.resize(12, 1);
  elastic_force.setZero();

  VectorXd line1;
  line1.resize(12, 1);
  line1.segment<3>(0) = tetrahedron.v_;
  for(unsigned int i = 0; i < 3; i++){
    line1.segment<3>((i+1)*3) = tetrahedron.p_inverse_.row(i).transpose();
  }

  elastic_force = 2*tetrahedron.lambda_*line1*(-3 + delta_0_);

  Matrix3d matrix1 = delta_;
  for(unsigned int i = 0; i < 3; i++){
    matrix1(i, i) *=2;
  }

  elastic_force.segment<3>(0) += tetrahedron.miu_*matrix1*tetrahedron.v_;
  for(unsigned int i = 0; i < 3; i++){
    elastic_force.segment<3>((i+1)*3) += tetrahedron.miu_*matrix1*tetrahedron.p_inverse_.row(i).transpose();
  }
  // cout << "===energy===" << endl;
  // cout << calculate_potential_energy(tetrahedron) << endl;
  // cout << "===force===" << endl;
  // cout << elastic_force << endl;
  // cout << "=== end ===" << endl;
  return elastic_force;
}

MatrixXd TetrahedronPotentialEnergyCalculator::calculate_delta_delta_potential_energy(Tetrahedron &tetrahedron){
  MatrixXd stiff;
  stiff.resize(12, 12);
  stiff.setZero();
  stiff.block(0, 0, 3, 3) = (2*tetrahedron.lambda_ + tetrahedron.miu_)*tetrahedron.v_*tetrahedron.v_.transpose()
    + (tetrahedron.miu_*tetrahedron.v_.transpose()*tetrahedron.v_)(0,0)*Matrix3d::Identity(); 
  for(unsigned int i = 1; i < 4; i++){
    stiff.block(0, i*3, 3, 3) = tetrahedron.miu_*tetrahedron.p_inverse_.row(i - 1).transpose()*tetrahedron.v_.transpose()
      + 2*tetrahedron.lambda_*tetrahedron.v_*tetrahedron.p_inverse_.row(i - 1)
      + (tetrahedron.miu_*tetrahedron.p_inverse_.row(i - 1)*tetrahedron.v_)(0, 0)*Matrix3d::Identity();
    //stiff.block(i*3, 0, 3, 3) = stiff.block(0, i*3, 3, 3).transpose();
    stiff.block(i*3, 0, 3, 3) = tetrahedron.miu_*tetrahedron.v_*tetrahedron.p_inverse_.row(i - 1)
      + 2*tetrahedron.lambda_*tetrahedron.p_inverse_.row(i - 1).transpose()*tetrahedron.v_.transpose()
      + (tetrahedron.miu_*tetrahedron.p_inverse_.row(i - 1)*tetrahedron.v_)(0, 0)*Matrix3d::Identity();
  }
  for(unsigned int i = 1; i < 4; i++){
    for(unsigned int j = 1; j < 4; j++){
      stiff.block(i*3, j*3, 3, 3) = 2*tetrahedron.lambda_*tetrahedron.p_inverse_.row(i - 1).transpose()*tetrahedron.p_inverse_.row(j - 1)
        + tetrahedron.miu_*tetrahedron.p_inverse_.row(j - 1).transpose()*tetrahedron.p_inverse_.row(i - 1) +
        + (tetrahedron.miu_*tetrahedron.p_inverse_.row(i - 1)*tetrahedron.p_inverse_.row(j - 1).transpose())(0,0)*Matrix3d::Identity();
    }
  }
  // cout << "miu : " << tetrahedron.miu_ << endl;
  // cout << "lambda : " << tetrahedron.lambda_ << endl;
  // cout << "p_inverse : " << tetrahedron.p_inverse_.transpose() << endl;
  // cout << "v_: " << tetrahedron.v_.transpose() << endl;

  return stiff;
}