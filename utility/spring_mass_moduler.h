#ifndef _SPRING_MASS_MODULER_H_
#define _SPRING_MASS_MODULER_H_

#include <Eigen/Dense>
#include <vector>
#include <utility>
#include "../model/spring.h"
#include "../model/point.h"


namespace zyclincoln{

  class SpringMassModuler{
  public:
    void AddPoint(std::vector<Eigen::Vector3d> &positions);
    void AddSpring(std::vector<std::pair<unsigned int, unsigned int>> pairs);
    void AddConstraint(const char *path);
    std::vector<Point> point();
    std::vector<Spring> spring();
  private:
    std::vector<unsigned int> index_of_static_points_;
    std::vector<Point> points_;
    std::vector<Spring> springs_;
  }

  inline std::vector<Point> SpringMassModuler::point(){
    return points_;
  }

  inline std::vector<Spring> SpringMassModuler::spring(){
    return springs_;
  }

}

#endif