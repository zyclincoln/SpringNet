#ifndef _POINT_H_
#define _POINT_H_

namespace zyclincoln{

  struct Point{
  public:
    Point(Eigen::Vector3d &position, Eigen::Vector3d &speed, double mass);
    Eigen::Vector3d position_;
    Eigen::Vector3d speed_;
    double mass_;
  }

}

#endif