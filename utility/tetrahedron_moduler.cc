#include "tetrahedron_moduler.h"
#include <fstream>
#include <iostream>
using namespace std;
using namespace Eigen;
using namespace zyclincoln;

void TetrahedronModuler::AddPoint(const vector<Vector3d> &positions){
  for(unsigned int i = 0; i < positions.size(); i++){
    points_.push_back(Point(positions[i], Vector3d(0, 0, 0), 10));
  }
}

void TetrahedronModuler::AddTetrahedron(const vector<vector<unsigned int>> &points_index){
  for(unsigned int i = 0; i < points_index.size(); i++){
    tetrahedrons_.push_back(Tetrahedron(0.2, 20, points_index[i]));
  }
}

void TetrahedronModuler::AddConstraint(const string &path){
  ifstream constraint_file;
  constraint_file.open(path);

  unsigned int static_point;
  unsigned int number;

  constraint_file >> number;

  for(unsigned int i = 0; i < number; i++){
    constraint_file >> static_point;
    index_of_static_points_.push_back(static_point - 1);
  }

  constraint_file.close();
}