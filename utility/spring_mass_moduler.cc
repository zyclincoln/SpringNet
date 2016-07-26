#include <assert.h>
#include <fstream>
#include "spring_mass_moduler.h"

using namespace std;
using namespace Eigen;
using namespace zyclincoln;

void SpringMassModuler::AddPoint(const vector<Vector3d> &positions){
  for(unsigned int i = 0; i < positions.size(); i++){
    points_.push_back(Point(positions[i], Vector3d(0, 0, 0), 10));
  }
}

void SpringMassModuler::AddSpring(const vector<pair<unsigned int, unsigned int>> &pairs, const vector<double> &edges_length){
  assert(pairs.size() == edges_length.size());

  for(unsigned int i = 0; i < pairs.size(); i++){
    springs_.push_back(Spring(40, edges_length[i], pairs[i]));
  }
}

void SpringMassModuler::AddConstraint(const string &path){
  ifstream constraint_file;
  constraint_file.open(path);

  unsigned int static_point;
  unsigned int number;

  constraint_file >> number;

  for(unsigned int i = 0; i < number; i++){
    constraint_file >> static_point;
    index_of_static_points_.push_back(static_point-1);  //index begin at 0
  }
  
  constraint_file.close();
}