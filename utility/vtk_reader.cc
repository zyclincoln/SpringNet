#include "vtk_reader.h"
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;
using namespace zyclincoln;
using namespace Eigen;

VTKReader::VTKReader(){

}

void VTKReader::ParseModel(){
  fstream fin;
  fin.open(file_name_);
  string line;

  while(1){
    fin>>line;
    if(fin.fail()){
      break;
    }
    if(line == string("POINTS")){
      unsigned int num;
      fin>>num;
      fin>>line;
      ParsePoint(fin, num);
    }
    else if(line == string("CELLS")){
      unsigned int num;
      fin>>num;
      fin>>line;
      ParseTetrahedron(fin, num);
    }
  }

  fin.close();

  cout<<"=======finish parse========"<<endl;
}

void VTKReader::ParsePoint(fstream &fin, unsigned int num){
  points_position_.clear();

  for(unsigned int i = 0; i < num; i++){
    double x, y, z;
    fin >> x >> y >> z;
    points_position_.push_back(Vector3d(x, y, z));
  
    cout<<"# point: "<<x<<" "<<y<<" "<<z<<endl;
  }

}

void VTKReader::ParseTetrahedron(fstream &fin, unsigned int num){
  tetrahedrons_point_.clear();

  for(unsigned int i = 0; i < num; i++){
    unsigned int id0, id1, id2, id3;
    fin>>id0>>id0>>id1>>id2>>id3;
    vector<unsigned int> index;
    index.push_back(id0);
    index.push_back(id1);
    index.push_back(id2);
    index.push_back(id3);

    tetrahedrons_point_.push_back(index);

    cout<<"# tetrahedron: "<<" "<<id0<<" "<<id1<<" "<<id2<<" "<<id3<<endl;
  }
}