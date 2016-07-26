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
    void AddPoint(const std::vector<Eigen::Vector3d> &positions);
    void AddSpring(const std::vector<std::pair<unsigned int, unsigned int>> &pairs, const std::vector<double> &edges_length);
    void AddConstraint(const std::string &path);
    std::vector<Point> point();
    std::vector<Spring> spring();
    std::vector<unsigned int> index_of_static_points();
  private:
    std::vector<unsigned int> index_of_static_points_;
    std::vector<Point> points_;
    std::vector<Spring> springs_;
  };

  inline std::vector<Point> SpringMassModuler::point(){
    return points_;
  }

  inline std::vector<Spring> SpringMassModuler::spring(){
    return springs_;
  }

  inline std::vector<unsigned int> SpringMassModuler::index_of_static_points(){
    return index_of_static_points_;
  }

}

#endif