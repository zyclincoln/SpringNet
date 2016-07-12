#include "obj_reader.h"
#include <iostream>
#include <sstream>

using namespace std;
using namespace Eigen;
using namespace zyclincoln;

using index = unsigned int;

OBJReader::OBJReader(const char* path):
	file_name_(path){

}

void OBJReader::ParseModel(){
	fstream fin;
	fin.open(file_name_);
	string line;
	while(getline(fin, line)){
		switch(line[0]){
			case 'v':{
				parse_point(line.substr(2));
				break;
			}
			case 'f':{
				parse_edge(line.substr(2));
				break;
			}
			default:{
				break;
			}
		}
	}
	fin.close();
}

void OBJReader::ParsePoint(string &line){
	double x,y,z;
	stringstream sstream;
	sstream << line;
	sstream >> x >> y >> z;
	points_position_.push_back(Vector3d(x, y, z));
}

void OBJReader::ParseEdge(string &line){
	index p0, p1, p2;
	stringstream sstream;

	sstream << line;
	sstream >> p0 >> p1 >> p2;
	p0--; p1--; p2--;	//index begins at 0

	if(unique_set_.find(pair<index, index>(p0, p1)) == unique_set_.end()) {
		unique_set_.insert(pair<index, index>(p0, p1));
		unique_set_.insert(pair<index, index>(p1, p0));
		unique_pair_.push_back(pair<index, index>(p0, p1));
		edges_length_.push_back((points_position_[p0]-points_position_[p1]).norm());
	}
	if(unique_set_.find(pair<index, index>(p0, p2)) == unique_set_.end()) {
		unique_set_.insert(pair<index, index>(p0, p2));
		unique_set_.insert(pair<index, index>(p2, p0));
		unique_pair_.push_back(pair<index, index>(p0, p2));
		edges_length_.push_back((points_position_[p0]-points_position_[p2]).norm());
	}
	if(unique_set_.find(pair<index, index>(p1, p2)) == unique_set_.end()) {
		unique_set_.insert(pair<index, index>(p1, p2));
		unique_set_.insert(pair<index, index>(p2, p1));
		unique_pair_.push_back(pair<index, index>(p1, p2));
		edges_length_.push_back((points_position_[p1]-points_position_[p2]).norm());
	}
}