#ifndef _VOL_2_VTK_H_
#define _VOL_2_VTK_H_

#include <string>
#include <iostream>
#include <fstream>
#include <vector>

namespace zyclincoln{
  void vol2vtk(std::string volpath, std::string vtkpath){
    std::ifstream fin;
    std::ofstream fout;
    fin.open(volpath);
  
    // scan file
    std::string st;
    std::vector<unsigned int> tetrahedron;
    std::vector<double> point;
    while(1){
      fin >> st;
      if(fin.fail()){
        break;
      }
      if(st == std::string("volumeelements")){
        unsigned int num;
        fin >> num;
        for(unsigned int i = 0; i < num; i++){
          unsigned int i0, i1, i2, i3;
          fin >> i0 >> i0;
          if(i0 != 4){
            std::cerr << "[ERROR] not tetrahedron!" << std::endl;
            fin.close();
            return;
          }
          fin >> i0 >> i1 >> i2 >> i3;
          tetrahedron.push_back(i0);
          tetrahedron.push_back(i1);
          tetrahedron.push_back(i2);
          tetrahedron.push_back(i3); 
        }
      }
  
      if(st == std::string("points")){
        unsigned int num;
        fin >> num;
        for(unsigned int i = 0; i < num; i++){
          double x, y, z;
          fin >> x >> y >> z;
          point.push_back(x);
          point.push_back(y);
          point.push_back(z); 
        }
      }
    }
    fin.close();
  
    // output result
    fout.open(vtkpath);
    fout << "# vtk DataFile Version 2.0\nTRI\nASCII\n\nDATASET UNSTRUCTURED_GRID\n";
    fout << "POINTS " << point.size()/3 << " float\n";
    for(unsigned int i = 0; i < point.size()/3; i++){
      fout << point[i*3] << " " << point[i*3+1] << " " << point[i*3+2] << "\n";
    }
  
    fout << "CELLS " << tetrahedron.size()/4 << " " << tetrahedron.size()/4*5 << std::endl;
    for(unsigned int i = 0; i < tetrahedron.size()/4; i++){
      fout << 4 << " " <<  tetrahedron[i*4+0]-1 << " " <<  tetrahedron[i*4+1]-1 <<  " " << tetrahedron[i*4+2]-1 << " " << tetrahedron[i*4+3]-1 << std::endl;
    }
  
    fout << "CELL_TYPES " << tetrahedron.size()/4 << std::endl;
    for(unsigned int i = 0; i < tetrahedron.size()/4; i++){
      fout << 10 << std::endl;
    }
  }
}

#endif