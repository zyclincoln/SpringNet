#ifndef _TETRAHEDRON_MODULER_H_
#define _TETRAHEDRON_MODULER_H_

#include <Eigen/Dense>
#include <vector>
#include <string>
#include "../model/point.h"
#include "../model/tetrahedron.h"

namespace zyclincoln{

  class TetrahedronModuler{
  public:
    void AddPoint(const std::vector<Eigen::Vector3d> &positions);
    void AddTetrahedron(const std::vector<std::vector<unsigned int>> &points_index);
    void AddConstraint(const std::string &path);
    std::vector<Point> point();
    std::vector<Tetrahedron> tetrahedron();
    std::vector<unsigned int> index_of_static_points();
  private:
    std::vector<unsigned int> index_of_static_points_;
    std::vector<Point> points_;
    std::vector<Tetrahedron> tetrahedrons_;
  };

  inline std::vector<unsigned int> TetrahedronModuler::index_of_static_points(){
    return index_of_static_points_;
  }

  inline std::vector<Point> TetrahedronModuler::point(){
    return points_;
  }

  inline std::vector<Tetrahedron> TetrahedronModuler::tetrahedron(){
    return tetrahedrons_;
  }

}


#endif