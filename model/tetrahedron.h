#ifndef _TETRAHEDRON_H_
#define _TETRAHEDRON_H_

#include <vector>
#include <Eigen/Dense>
#include "point.h"

namespace zyclincoln{

  struct Tetrahedron{
  public:
    Tetrahedron(const double poison, const double young, const std::vector<unsigned int> points_index);
    PreCompute(const std::vector<Point> &points);
    ComputeDx(const std::vector<Point> &points);
    double poison_;
    double young_;
    std::vector<unsigned int> points_index_;
    Eigen::MatrixXd dx_;
    Eigen::MatrixXd p_inverse_;
    Eigen::VectorXd v_;
    double lambda = 0;
    double miu = 0;
  };

}

#endif