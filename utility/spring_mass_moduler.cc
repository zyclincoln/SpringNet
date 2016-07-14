#include <assert.h>
#include <fstream>
#include "spring_mass_moduler.h"

using namespace std;
using namespace Eigen;
using namespace zyclincoln;

using index = unsigned int;

void SpringMassModuler::AddPoint(vector<Vector3d> &positions){
  for(index i = 0; i < positions.size(); i++){
    points_.push_back(Point(positions[i], Vector3d(0, 0, 0), 1));
  }
}

void SpringMassModuler::AddSpring(vector<pair<index, index>> pairs, vector<double> edges_length){
  assert(pairs.size() == edges_length.size());

  for(index i = 0; i < pairs.size(); i++){
    springs_.push_back(Spring(1, edges_length[i], pairs[i]));
  }
}

void SpringMassModuler::AddConstraint(const char *path){
  ifstream constraint_file;
  constraint_file.open(path);

  index static_point;
  unsigned int number;

  constraint_file >> number;

  for(index i = 0; i < number; i++){
    constraint_file >> static_point;
    index_of_static_points_.push_back(static_point-1);  //index begin at 0
  }
  
  constraint_file.close();
}