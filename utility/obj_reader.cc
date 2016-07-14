#include "obj_reader.h"
#include <iostream>
#include <sstream>

using namespace std;
using namespace Eigen;
using namespace zyclincoln;

OBJReader::OBJReader(){

}

void OBJReader::ParseModel(){
	fstream fin;
	fin.open(file_name_);
	string line;
	while(getline(fin, line)){
		switch(line[0]){
			case 'v':{
				ParsePoint(line.substr(2));
				break;
			}
			case 'f':{
				ParseEdge(line.substr(2));
				break;
			}
			default:{
				break;
			}
		}
	}
	fin.close();
}

void OBJReader::ParsePoint(string line){
	double x,y,z;
	stringstream sstream;
	sstream << line;
	sstream >> x >> y >> z;
	points_position_.push_back(Vector3d(x, y, z));
}

void OBJReader::ParseEdge(string line){
	unsigned int p0, p1, p2;
	stringstream sstream;

	sstream << line;
	sstream >> p0 >> p1 >> p2;
	p0--; p1--; p2--;	//index begins at 0

	if(unique_set_.find(pair<unsigned int, unsigned int>(p0, p1)) == unique_set_.end()) {
		unique_set_.insert(pair<unsigned int, unsigned int>(p0, p1));
		unique_set_.insert(pair<unsigned int, unsigned int>(p1, p0));
		unique_pair_.push_back(pair<unsigned int, unsigned int>(p0, p1));
		edges_length_.push_back((points_position_[p0]-points_position_[p1]).norm());
	}
	if(unique_set_.find(pair<unsigned int, unsigned int>(p0, p2)) == unique_set_.end()) {
		unique_set_.insert(pair<unsigned int, unsigned int>(p0, p2));
		unique_set_.insert(pair<unsigned int, unsigned int>(p2, p0));
		unique_pair_.push_back(pair<unsigned int, unsigned int>(p0, p2));
		edges_length_.push_back((points_position_[p0]-points_position_[p2]).norm());
	}
	if(unique_set_.find(pair<unsigned int, unsigned int>(p1, p2)) == unique_set_.end()) {
		unique_set_.insert(pair<unsigned int, unsigned int>(p1, p2));
		unique_set_.insert(pair<unsigned int, unsigned int>(p2, p1));
		unique_pair_.push_back(pair<unsigned int, unsigned int>(p1, p2));
		edges_length_.push_back((points_position_[p1]-points_position_[p2]).norm());
	}
}