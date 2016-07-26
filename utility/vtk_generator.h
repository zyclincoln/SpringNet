#ifndef _VTK_GENERATOR_H_
#define _VTK_GENERATOR_H_

#include <string>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include "../model/point.h"
#include "../model/tetrahedron.h"
#include "../model/fem_system.h"

using std::string;
using std::ofstream;
using std::vector;
using std::endl;

namespace zyclincoln{
  void linear_fem_system_to_vtk(LinearFEMSystem &system, string &path){
    ofstream fout;
    fout.open(path);
    fout << "# vtk DataFile Version 2.0\nTRI\nASCII\n\nDATASET UNSTRUCTURED_GRID\n";
  
    vector<Point> point_set = system.points();
    unsigned int point_size = point_set.size();
    fout << "POINTS " << point_size << " float\n";
    for(unsigned int i = 0; i < point_size; i++){
      fout << point_set[i].position_(0, 0) <<" " << point_set[i].position_(1, 0) << " " << point_set[i].position_(2, 0)<< "\n";
    }
  
    vector<Tetrahedron> tetrahedron_set = system.tetrahedrons();
    unsigned int tetrahedron_size = tetrahedron_set.size();
    fout << "CELLS " << tetrahedron_size << " " << tetrahedron_size*5 << endl;
    for(unsigned int i = 0; i < tetrahedron_size; i++){
      fout << 4 << " " <<  tetrahedron_set[i].points_index_[0] << " " <<  tetrahedron_set[i].points_index_[1]
        <<  " " << tetrahedron_set[i].points_index_[2] << " " << tetrahedron_set[i].points_index_[3] << endl;
    }
  
    fout << "CELL_TYPES " << tetrahedron_size << endl;
    for(unsigned int i = 0; i <tetrahedron_size; i++){
      fout << 10 << endl;
    }
  }
}
#endif