#ifndef _VTK_READER_H_
#define _VTK_READER_H_

#include <vector>
#include <Eigen/Dense>
#include <fstream>
#include <string>

namespace zyclincoln{

  class VTKReader{
  public:
    VTKReader();
    void set_file_name(const std::string path);
    std::vector<Eigen::Vector3d>& points_position();
    unsigned int points_num();
    std::vector<std::vector<unsigned int>>& tetrahedrons_point();
    unsigned int tetrahedrons_num();
    void ParseModel();
  private:
    void ParsePoint(std::fstream &fin, unsigned int num);
    void ParseTetrahedron(std::fstream &fin, unsigned int num);
    std::vector<Eigen::Vector3d> points_position_;
    std::vector<std::vector<unsigned int>> tetrahedrons_point_;
    std::string file_name_;
  };

  inline std::vector<Eigen::Vector3d>& VTKReader::points_position(){
    return points_position_;
  }

  inline void VTKReader::set_file_name(const std::string path){
    file_name_ = path;
  }

  inline std::vector<std::vector<unsigned int>>& VTKReader::tetrahedrons_point(){
    return tetrahedrons_point_;
  }

  inline unsigned int VTKReader::points_num(){
    return points_position_.size();
  }

  inline unsigned int VTKReader::tetrahedrons_num(){
    return tetrahedrons_point_.size();
  }

}

#endif