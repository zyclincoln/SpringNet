#ifndef _NET_READER_H_
#define _NET_READER_H_

#include <vector>
#include <string>
#include <fstream>
#include <Eigen/Dense>

class Net_Reader{
public:
	Net_Reader(const char* path);
	unsigned int length();
	unsigned int width();
	std::vector<double> mass();
	std::vector<Eigen::Vector3d> points_position();
	std::vector<Eigen::Vector3d> points_speed();
	std::vector<unsigned int> index_of_static_points();
	void read_file();
private:
	void read_size();
	void read_points();
	unsigned int length_;
	unsigned int width_;
	std::vector<double> mass_;
	std::vector<Eigen::Vector3d> points_position_;
	std::vector<Eigen::Vector3d> points_speed_;
	std::vector<unsigned int> index_of_static_points_;
	std::ifstream fin;
};

inline unsigned int Net_Reader::length(){
	return length_;
}

inline unsigned int Net_Reader::width(){
	return width_;
}

inline std::vector<double> Net_Reader::mass(){
	return mass_;
}

inline std::vector<Eigen::Vector3d> Net_Reader::points_position(){
	return points_position_;
}

inline std::vector<Eigen::Vector3d> Net_Reader::points_speed(){
	return points_speed_;
}

inline std::vector<unsigned int> Net_Reader::index_of_static_points(){
	return index_of_static_points_;
}

#endif