#include "utility/vol2vtk.h"
#include <string>
#include <iostream>

using namespace std;
using namespace zyclincoln;

int main(int argc, char** argv){
  if(argc != 3){
    cerr << "[ERROR] parameter is not enough!" << endl;
    return -1;
  }
  string vol_path = string(argv[1]);
  string vtk_path = string(argv[2]);
  vol2vtk(vol_path, vtk_path);
  return 0;
}
