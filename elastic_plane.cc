#include <assert.h>
#include <math.h>
#include "Elastic_plane.h"

using namespace Eigen;
using namespace std;

void ElasticPlane::default_spring_layout(){
	unsigned int index = 0;
	
	for(unsigned int i = 0; i <= cwidth_; i++){
		for(unsigned int j = 0; j < clength_; j++){
			springs_[index].stiff_ = 1;
			springs_[index].length_ = 1;
			springs_[index].between_.first = i*clength_ + j;
			springs_[index].between_.second = i*clength_ + j + 1;
			index++;
		}
	}

	for(unsigned int i = 0; i <= clength_; i++){
		for(unsigned int j = 0; j < cwidth_; j++){
			springs_[index].stiff_ = 1;
			springs_[index].length_ = 1;
			springs_[index].between_.first = j*clength_ + i;
			springs_[index].between_.second = (j+1)*clength_ + i;
			index++; 
		}
	}

	for(unsigned int i=0; i < ciwth_ ;i++){
		for(unsigned int j = 0; j < clength_; j++){
			springs_[index].stiff_ = 1;
			springs_[index].length_ = 1;
			springs_[index].between_.first = i*clength_ + j;
			springs_[index].between_.second = (i+1)*clenght_ + j + 1;
		}
	}
}

void ElasticPlane::setup_matrix(){
	mass_matrix_.resize((clength_+1)*(cwidth_+1)*3, (clength_+1)*(cwidth_+1)*3);
	stiffness_matrix_.resize((clength_+1)*(cwidth_+1)*3, (clength_+1)*(cwidth_+1)*3);
	points_position_.resize((clength_+1)*(cwidth_+1)*3, 1);
	points_speed_.resize((clength_+1)*(cwidth_+1)*3, 1);
	springs_.resize(springs_);

	stiffness_matrix_.setZero();
	mass_matrix_.setZero();
	points_speed_.setZero();
	points_position_.setZero();
}

void ElasticPlane::next_frame(){
	stiffness_matrix_.setZero();
	VectorXd force = VectorXd::Zero((clength_+1)*(cwidth_+1)*3, 1);

	for(unsigned int i = 0; i < springs_.size(); i++){
		unsigned index0 = springs_[i].between.first;
		unsigned index1 = springs_[i].between.second;
		Vector3d p0 = points_position_.block(index0*3, 0, 3, 1);
		Vector3d p1 = points_position_.block(index1*3, 0, 3, 1);

		Matrix3d sub_stiff = Matrix3d::Identity()*(1 - 1/springs_[i].length);
		sub_stiff += (p1 - p0)*(p1 - p0).transpose()/pow(springs_[i].length,3);

		for(int i = 0; i < 3; i++){
			stiffness_matrix_.block(index0*3, index0*3, 3, 3) += sub_stiff;
			stiffness_matrix_.block(index1*3, index1*3, 3, 3) += sub_stiff;
			stiffness_matrix_.block(index0*3, index1*3, 3, 3) -= sub_stiff;
			stiffness_matrix_.block(index1*3, index0*3, 3, 3) -= sub_stiff;
		}

		force(index0*3, 0, 3, 1) += ((p1 - p0).norm() - springs_[i].length)*springs_[i].stiff_*(p1 - p0).normalized();
		force(index1*3, 0, 3, 1) += ((p0 - p1).norm() - springs_[i].length)*springs_[i].stiff_*(p0 - p1).normalized();
	}

	MatrixXd A = MatrixXd::Zero((clength_+1)*(cwidth_+1)*3);
	A = mass_matrix_ + pow(time_step_ms_, 2)*stiffness_matrix_;
	VectorXd b = VectorXd::Zero((clength_+1)*(cwidth_+1)*3, 1);
	b = time_step_ms_*(force - stiffness_matrix_*time_step_ms_*points_speed_);

	set<unsigned int> index_need_update;
	for(unsigned int i=0; i < (clength_+1)*(cwidth_+1); i++){
		if(static_points_.find(i)!=static_points_.end()){
			set.insert(i*3);
			set.insert(i*3+1);
			set.insert(i*3+2);
		}
	}

	MatrixXd shrinked_A;
	VectorXd shrinked_b;

	shrinked_A = shrink_matrix(A, index_need_update, index_need_update);
	shrinked_b = shrink_colvector(b, index_need_update);

	VectorXd shrinked_delta_v;
	shrinked_delta_v = shrinked_A.colPivHouseholderQr().solve(shrinked_b);

	VectorXd delta_v = VectorXd::Zero((clength_+1)*(cwidth_+1)*3, 1);
	map_to_original_colvector(shrinked_delta_v, index_need_update, delta_v);

	points_speed_ += delta_v;
	points_position_ += delta_v*time_step_ms_;
}

void ElasticPlane::add_static_points(std::vector<unsigned int> index_of_static_points){
	for(unsigned int i = 0; i < index_of_static_points.size(); i++){
		static_points_.insert(index_of_static_points[i]);
	}
}

void ElasticPlane::add_points(const vector<double> &mass, const vector<Vertex3d> &position, const vector<Vertex3d> &speed){
	assert(mass.size() == position.size());
	assert(mass.size() == (clength_+1)*(cwidth_+1));

	for(unsigned int i = 0; i < mass.size(); i++){
		mass_.block(i*3, i*3, 3, 3) = Matrix3d::Identity()*mass[i];
		points_position_.(i*3, 1, 3, 1) = position[i];
		points_speed_.(i*3, 1, 3, 1) = speed[i];
	}
}

void ElasticPlane::add_point(const unsigned int index, const double mass, const Vector3d &positon, const Vector3d &speed){
	assert(index < (clength_+1)*(cwidth_+1));
	
	mass_.block(index*3, index*3, 3, 3) = Matrix3d::Identity()*mass;
	points_position_.(i*3, 1, 3, 1) = position;
	points_speed_.(i*3, 1, 3, 1) = speed;
}

void ElasticPlane::add_springs(const vector<double> &stiffness, const vector<double> &length, const vector<pair<int, int> > &between){
	assert(stiffness.size() == length.size());
	assert(stiffness.size() == betweeen.size());
	assert(stiffness.size() == springs_.size());

	for(unsigned int i = 0; i < stiffness.size(); i++){
		springs_[i].stiff_ = stiffness[i];
		springs_[i].length_ = length[i];
		springs_[i].between_ = between[i];
	}
}

void ElasticPlane::add_spring(const unsigned int index, const double stiffness, const double length, const pair<int, int> &between){
	assert(index < springs_.size());

	springs_[index].stiff_ = stiffness;
	springs_[index].length_ = length;
	springs_[index].between_ = between;
}