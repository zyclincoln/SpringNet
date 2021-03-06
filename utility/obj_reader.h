#ifndef _OBJ_READER_H_
#define _OBJ_READER_H_

#include <vector>
#include <utility>
#include <set>
#include <string>
#include <fstream>
#include <Eigen/Dense>

namespace zyclincoln{

	class OBJReader{
	public:
		OBJReader();
		void set_file_name(const std::string path);
		std::vector<Eigen::Vector3d>& points_position();
		std::vector<std::pair<unsigned int, unsigned int> > unique_pairs();
		unsigned int points_num();
		unsigned int unique_pairs_num();
		std::vector<double> edges_length();
		void ParseModel();
	private:
		void ParsePoint(std::string line);
		void ParseEdge(std::string line);
		std::vector<Eigen::Vector3d> points_position_;
		std::vector<double> edges_length_;
		std::set<std::pair<unsigned int, unsigned int> > unique_set_;
		std::vector<std::pair<unsigned int, unsigned int> > unique_pair_;
		std::string file_name_;
	};

	inline std::vector<Eigen::Vector3d>& OBJReader::points_position(){
		return points_position_;
	}

	inline void OBJReader::set_file_name(const std::string path){
		file_name_ = path;
	}

	inline std::vector<std::pair<unsigned int, unsigned int> > OBJReader::unique_pairs(){
		return unique_pair_;
	}

	inline unsigned int OBJReader::points_num(){
		return points_position_.size();
	}

	inline unsigned int OBJReader::unique_pairs_num(){
		return edges_length_.size();
	}

	inline std::vector<double> OBJReader::edges_length(){
		return edges_length_;
	}

}

#endif