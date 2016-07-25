#include "../utility/vtk_reader.h"

#include <string>
using namespace std;
using namespace zyclincoln;

int main(){
  VTKReader vtk_reader;
  vtk_reader.set_file_name(string("../build/beam.vtk"));
  vtk_reader.ParseModel();
  return 0;
}