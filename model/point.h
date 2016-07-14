#ifndef _POINT_H_
#define _POINT_H_

#include <Eigen/Dense>

namespace zyclincoln{

  struct Point{
  public:
    Point(const Eigen::Vector3d &position, const Eigen::Vector3d &speed, const double mass);
    Eigen::Vector3d position_;
    Eigen::Vector3d speed_;
    double mass_;
  };

}

#endif