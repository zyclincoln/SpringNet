#include <iostream>
#include "net_reader.h"

using namespace std;

Net_Reader::Net_Reader(const char* path){
	fin.open(path);
}

void Net_Reader::read_file(){
	read_size();
	read_points();
	fin.close();
}

void Net_Reader::read_size(){
	fin>>length_>>width_;	
}

void Net_Reader::read_points(){
	double mass;
	Eigen::Vector3d position, speed;
	unsigned int static_index;

	for(unsigned int i = 0; i < (length_+1)*(width_+1); i++){
		fin >> mass >> position(0,0) >> position(1,0) >> position(2,0) >> speed(0,0) >> speed(1,0) >> speed(2,0);
		mass_.push_back(mass);
		points_position_.push_back(position);
		points_speed_.push_back(speed);
	}

	unsigned int static_point_num;
	fin>>static_point_num;
	for(unsigned int i=0; i<static_point_num; i++){
		fin >> static_index;
		index_of_static_points_.push_back(static_index);
	}
}