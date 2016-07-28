#ifndef _TETRAHEDRON_H_
#define _TETRAHEDRON_H_

#include <vector>
#include <Eigen/Dense>
#include "point.h"

namespace zyclincoln{

  struct Tetrahedron{
  public:
    Tetrahedron(const double poison, const double young, const std::vector<unsigned int> &points_index);
    void PreCompute(const std::vector<Point> &points);
    void ComputeDx(const std::vector<Point> &points);
    void ComputeR();    // make sure the Dx has been updated
    double poison_;
    double young_;
    std::vector<unsigned int> points_index_;
    Eigen::MatrixXd dx_;
    Eigen::MatrixXd p_inverse_;
    Eigen::VectorXd v_;
    double lambda_ = 0;
    double miu_ = 0;
    Eigen::MatrixXd r_;
  };

}

#endif