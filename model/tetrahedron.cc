#include "tetrahedron.h"
#include <assert.h>
#include <iostream>
using namespace std;
using namespace zyclincoln;
using namespace Eigen;

Tetrahedron::Tetrahedron(const double poison, const double young, const vector<unsigned int> &points_index):
  poison_(poison),
  young_(young){

  assert(points_index.size() == 4);

  points_index_ = points_index;
  p_inverse_ = MatrixXd::Zero(3, 3);
  v_ = Vector3d::Zero();

  miu_ = young_/(2*(1+poison_));
  lambda_ = poison_*young_/((1 + poison_)*(1 - 2*poison_));
}

void Tetrahedron::PreCompute(const std::vector<Point> &points){
  assert(points.size() == 4);

  MatrixXd position = MatrixXd::Zero(3, 3);

  for(unsigned int i = 0; i < 3; i++){
    position.col(i) = points[i + 1].position_ - points[0].position_;
  }

  p_inverse_ = position.inverse();

  v_ = -p_inverse_.row(0).transpose() - p_inverse_.row(1).transpose() - p_inverse_.row(2).transpose();

}

void Tetrahedron::ComputeDx(const std::vector<Point> &points){
  assert(points.size() == 4);

  dx_ = MatrixXd::Zero(3, 3);

  for(unsigned int i = 0; i < 3; i++){
    dx_.col(i) = points[i+1].position_ - points[0].position_;
  }
}