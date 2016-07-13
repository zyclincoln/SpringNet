#include "obj_reader.h"
#include <iostream>
#include <sstream>

using namespace std;
using namespace Eigen;

Obj_Reader::Obj_Reader(const char* path){
	fin.open(path);
}

void Obj_Reader::sparse_model(){
	string line;
	while(getline(fin, line)){
		switch(line[0]){
			case 'v':
				sparse_point(line.substr(2));
				break;
			case 'f':
				sparse_edge(line.substr(2));
				break;
			default:
				break;
		}
	}
	fin.close();

	// for(int i=0; i<mass_.size(); i++){
	// 	cout<<i<<" : mass "<<mass_[i]<<" position "<<points_position_[i].transpose()<<" speed "<<points_speed_[i].transpose()<<endl;
	// }

	// for(int i=0; i<springs_stiff_.size(); i++){
	// 	cout<<i<<" : stiff "<<springs_stiff_[i]<<" length "<<springs_length_[i]<<" between "<<springs_pair_[i].first<<" "<<springs_pair_[i].second<<endl;
	// }
}

void Obj_Reader::sparse_constraint(const char *path){
	ifstream constraint_file;
	constraint_file.open(path);

	unsigned int static_point;
	unsigned int number;

	constraint_file >> number;

	for(int i=0; i < number; i++){
		constraint_file >> static_point;
		index_of_static_points_.push_back(static_point-1);
	}
	
	constraint_file.close();
}

void Obj_Reader::sparse_point(string line){
	double x,y,z;
	stringstream sstream;
	sstream << line;
	sstream >> x >> y >> z;
	points_position_.push_back(Vector3d(x, y, z));
	points_speed_.push_back(Vector3d(0, 0, 0));
	mass_.push_back(10);
}

void Obj_Reader::sparse_edge(string line){
	unsigned p0, p1, p2;
	stringstream sstream;

	sstream << line;
	sstream >> p0 >> p1 >> p2;
	p0--; p1--; p2--;
	if(springs_set_.find(pair<unsigned int, unsigned int>(p0, p1)) == springs_set_.end()) {
		springs_set_.insert(pair<unsigned int, unsigned int>(p0, p1));
		springs_set_.insert(pair<unsigned int, unsigned int>(p1, p0));
		springs_pair_.push_back(pair<unsigned int, unsigned int>(p0, p1));
		springs_stiff_.push_back(1);
		springs_length_.push_back((points_position_[p0]-points_position_[p1]).norm());
	}
	if(springs_set_.find(pair<unsigned int, unsigned int>(p0, p2)) == springs_set_.end()) {
		springs_set_.insert(pair<unsigned int, unsigned int>(p0, p2));
		springs_set_.insert(pair<unsigned int, unsigned int>(p2, p0));
		springs_pair_.push_back(pair<unsigned int, unsigned int>(p0, p2));
		springs_stiff_.push_back(1);
		springs_length_.push_back((points_position_[p0]-points_position_[p2]).norm());
	}
	if(springs_set_.find(pair<unsigned int, unsigned int>(p1, p2)) == springs_set_.end()) {
		springs_set_.insert(pair<unsigned int, unsigned int>(p1, p2));
		springs_set_.insert(pair<unsigned int, unsigned int>(p2, p1));
		springs_pair_.push_back(pair<unsigned int, unsigned int>(p1, p2));
		springs_stiff_.push_back(1);
		springs_length_.push_back((points_position_[p1]-points_position_[p2]).norm());
	}
}