#ifndef _OBJ_READER_H_
#define _OBJ_READER_H_

#include <vector>
#include <utility>
#include <set>
#include <string>
#include <fstream>
#include <Eigen/Dense>

class Obj_Reader{
public:
	Obj_Reader(const char* path);
	void sparse_model();
	void sparse_constraint(const char *path);
	std::vector<Eigen::Vector3d> points_position();
	std::vector<Eigen::Vector3d> points_speed();
	std::vector<double> mass();
	std::vector<unsigned int> index_of_static_points();
	std::vector<std::pair<unsigned int, unsigned int> > springs_pair();
	std::vector<double> springs_stiff();
	std::vector<double> springs_length();
	unsigned int points_num();
	unsigned int springs_num();
private:
	void sparse_point(std::string line);
	void sparse_edge(std::string line);
	std::vector<Eigen::Vector3d> points_position_;
	std::vector<Eigen::Vector3d> points_speed_;
	std::vector<double> mass_;
	std::vector<unsigned int> index_of_static_points_;
	std::set<std::pair<unsigned int, unsigned int> > springs_set_;
	std::vector<std::pair<unsigned int, unsigned int> > springs_pair_;
	std::vector<double> springs_stiff_;
	std::vector<double> springs_length_;
	std::ifstream fin;
};

inline std::vector<Eigen::Vector3d> Obj_Reader::points_position(){
	return points_position_;
}

inline std::vector<Eigen::Vector3d> Obj_Reader::points_speed(){
	return points_speed_;
}

inline std::vector<double> Obj_Reader::mass(){
	return mass_;
}

inline std::vector<unsigned int> Obj_Reader::index_of_static_points(){
	return index_of_static_points_;
}

inline std::vector<std::pair<unsigned int, unsigned int> > Obj_Reader::springs_pair(){
	return springs_pair_;
}

inline std::vector<double> Obj_Reader::springs_stiff(){
	return springs_stiff_;
}

inline std::vector<double> Obj_Reader::springs_length(){
	return springs_length_;
}

inline unsigned int Obj_Reader::points_num(){
	return mass_.size();
}

inline unsigned int Obj_Reader::springs_num(){
	return springs_length_.size();
}

#endif