#include <assert.h>
#include <math.h>
#include <iostream>
#include "elastic_plane.h"
#include "matrix_op.h"

using namespace Eigen;
using namespace std;

void ElasticPlane::generate_draw_line_index(){
	draw_line_index_ = VectorXd::Zero(springs_.size()*2, 1);

	for(unsigned int i = 0; i < springs_.size(); i++){
		draw_line_index_(i*2 + 0, 0) = springs_[i].between_.first;
		draw_line_index_(i*2 + 1, 0) = springs_[i].between_.second;
	}
}

void ElasticPlane::next_frame(){
	set<unsigned int> index_need_remove;
	for(unsigned int i=0; i < cpoints_; i++){
		if(static_points_.find(i)!=static_points_.end()){
			index_need_remove.insert(i*3);
			index_need_remove.insert(i*3+1);
			index_need_remove.insert(i*3+2);
		}
	}


	VectorXd pt = points_position_;

	for(int iter = 0; iter < 1 ; iter++) {
		stiffness_matrix_.setZero();
		VectorXd force = VectorXd::Zero(cpoints_*3, 1);
		for(unsigned int i = 0; i < springs_.size(); i++){
			unsigned index0 = springs_[i].between_.first;
			unsigned index1 = springs_[i].between_.second;
			Vector3d p0 = pt.segment<3>(index0*3);
			Vector3d p1 = pt.segment<3>(index1*3);

			Matrix3d sub_stiff = Matrix3d::Identity()*(1 - 1/springs_[i].length_);
			sub_stiff += (p1 - p0)*(p1 - p0).transpose()/pow(springs_[i].length_,3);
			stiffness_matrix_.block(index0*3, index0*3, 3, 3) += sub_stiff;
			stiffness_matrix_.block(index1*3, index1*3, 3, 3) += sub_stiff;
			stiffness_matrix_.block(index0*3, index1*3, 3, 3) -= sub_stiff;
			stiffness_matrix_.block(index1*3, index0*3, 3, 3) -= sub_stiff;
			force.segment<3>(index0*3) += ((p1 - p0).norm() - springs_[i].length_)*springs_[i].stiff_*(p1 - p0).normalized();
			force.segment<3>(index1*3) += ((p1 - p0).norm() - springs_[i].length_)*springs_[i].stiff_*(p0 - p1).normalized();
	
			// cout<<"===part force===\n"<<
			// 	"p0 : "<<p0.transpose()<<endl<<
			// 	"p1 : "<<p1.transpose()<<endl<<
			// 	"offset : "<<(p1-p0).norm()-springs_[i].length_<<endl<<
			// 	"direction : "<<(p1 - p0).normalized().transpose()<<endl<<
			// 	"force : "<<((p1 - p0).norm() - springs_[i].length_)*springs_[i].stiff_*(p1 - p0).normalized().transpose()<<endl;
		}
	
		for(unsigned int i=0; i < cpoints_; i++){
			force.segment<3>(i*3) += Vector3d(0, -9.8*0.001, 0);
		}
	
		// cout<<"===total force===\n"<<force.transpose()<<endl;
	
		MatrixXd A = MatrixXd::Zero(cpoints_*3, cpoints_*3);
		A = mass_matrix_ + pow(time_step_ms_, 2)*stiffness_matrix_;
		VectorXd b = VectorXd::Zero(cpoints_*3, 1);
		b = -mass_matrix_*(pt - points_position_ - time_step_ms_*points_speed_) + pow(time_step_ms_, 2) * force;
		// b = time_step_ms_*(force - stiffness_matrix_*time_step_ms_*points_speed_);
		MatrixXd shrinked_A;
		VectorXd shrinked_b;
	
		shrink_matrix(A, index_need_remove, index_need_remove, shrinked_A);
		shrink_colvector(b, index_need_remove, shrinked_b);
		// cout<<"===A===\n"<<A<<endl;
		// cout<<"===shrinked_A===\n"<<shrinked_A<<endl;
		// cout<<"===b===\n"<<b.transpose()<<endl;
		// cout<<"===shrinked_b===\n"<<shrinked_b.transpose()<<endl;
		VectorXd shrinked_delta_p = shrinked_A.colPivHouseholderQr().solve(shrinked_b);
		VectorXd delta_p = VectorXd::Zero(cpoints_*3, 1);
		map_to_original_colvector(shrinked_delta_p, index_need_remove, delta_p);
		// cout<<"===delta_p===\n"<<delta_p.transpose()<<endl;
		pt += delta_p;
		// cout<<"===pt===\n"<<pt.transpose()<<endl;
		//getchar();
	}

	points_speed_ = (pt - points_position_)/time_step_ms_;
	points_position_ = pt;
	// cout<<"===point position===\n"<<points_position_<<endl;
}

void ElasticPlane::next_frame_2D(){
	set<unsigned int> index_need_remove;
	for(unsigned int i=0; i < cpoints_; i++){
		if(static_points_.find(i)!=static_points_.end()){
			index_need_remove.insert(i*2);
			index_need_remove.insert(i*2+1);
		}
	}

	VectorXd pt = points_position_2D_;
	for(int iter = 0; iter < 1 ; iter++) {
		stiffness_matrix_2D_.setZero();
		VectorXd force = VectorXd::Zero(cpoints_*2, 1);
		for(unsigned int i = 0; i < springs_.size(); i++){
			unsigned index0 = springs_[i].between_.first;
			unsigned index1 = springs_[i].between_.second;
			Vector2d p0 = pt.segment<2>(index0*2);
			Vector2d p1 = pt.segment<2>(index1*2);


			Matrix2d sub_stiff = Matrix2d::Identity()*(1 - 1/springs_[i].length_);
			sub_stiff += (p1 - p0)*(p1 - p0).transpose()/pow(springs_[i].length_,3);
			stiffness_matrix_.block(index0*2, index0*2, 2, 2) += sub_stiff;
			stiffness_matrix_.block(index1*2, index1*2, 2, 2) += sub_stiff;
			stiffness_matrix_.block(index0*2, index1*2, 2, 2) -= sub_stiff;
			stiffness_matrix_.block(index1*2, index0*2, 2, 2) -= sub_stiff;
			force.segment<2>(index0*2) += ((p1 - p0).norm() - springs_[i].length_)*springs_[i].stiff_*(p1 - p0).normalized();
			force.segment<2>(index1*2) += ((p1 - p0).norm() - springs_[i].length_)*springs_[i].stiff_*(p0 - p1).normalized();
		}

	
		for(unsigned int i=0; i < cpoints_; i++){
			force.segment<2>(i*2) += Vector2d(0, -9.8*0.001);
		}

	
		MatrixXd A = MatrixXd::Zero(cpoints_*2, cpoints_*2);
		A = mass_matrix_2D_ + pow(time_step_ms_, 2)*stiffness_matrix_2D_;
		VectorXd b = VectorXd::Zero(cpoints_*2, 1);
		b = -mass_matrix_2D_*(pt - points_position_2D_ - time_step_ms_*points_speed_2D_) + pow(time_step_ms_, 2) * force;

		MatrixXd shrinked_A;
		VectorXd shrinked_b;
    
		shrink_matrix(A, index_need_remove, index_need_remove, shrinked_A);
		shrink_colvector(b, index_need_remove, shrinked_b);
		VectorXd shrinked_delta_p = shrinked_A.colPivHouseholderQr().solve(shrinked_b);
		VectorXd delta_p = VectorXd::Zero(cpoints_*2, 1);
		map_to_original_colvector(shrinked_delta_p, index_need_remove, delta_p);
		pt += delta_p;
	}

	points_speed_2D_ = (pt - points_position_2D_)/time_step_ms_;
	points_position_2D_ = pt;

	// cout<<"===position===\n"<<points_position_2D_<<endl;
}

void ElasticPlane::add_static_points(std::vector<unsigned int> index_of_static_points){
	for(unsigned int i = 0; i < index_of_static_points.size(); i++){
		static_points_.insert(index_of_static_points[i]);
	}
}

void ElasticPlane::add_points(const vector<double> &mass, const vector<Vector3d> &position, const vector<Vector3d> &speed){

	for(unsigned int i = 0; i < mass.size(); i++){
		mass_matrix_.block(i*3, i*3, 3, 3) = Matrix3d::Identity()*mass[i];
		points_position_.block(i*3, 0, 3, 1) = position[i];
		points_speed_.block(i*3, 0, 3, 1) = speed[i];
	}
}

void ElasticPlane::add_points_2D(const vector<double> &mass, const vector<Vector3d> &position, const vector<Vector3d> &speed){

	for(unsigned int i = 0; i < mass.size(); i++){
		mass_matrix_2D_.block(i*2, i*2, 2, 2) = Matrix2d::Identity()*mass[i];
		points_position_2D_.segment<2>(i*2) = position[i].segment<2>(0);
		points_speed_2D_.segment<2>(i*2) = speed[i].segment<2>(0);
	}
}

void ElasticPlane::add_point(const unsigned int index, const double mass, const Vector3d &position, const Vector3d &speed){
	assert(index < cpoints_);
	
	mass_matrix_.block(index*3, index*3, 3, 3) = Matrix3d::Identity()*mass;
	points_position_.block(index*3, 0, 3, 1) = position;
	points_speed_.block(index*3, 0, 3, 1) = speed;
}

void ElasticPlane::add_springs(const vector<double> &stiffness, const vector<double> &length, const vector<pair<unsigned int, unsigned int> > &between){
	assert(stiffness.size() == length.size());
	assert(stiffness.size() == between.size());

	for(unsigned int i = 0; i < stiffness.size(); i++){
		springs_.push_back(Spring(stiffness[i]*200, length[i], between[i]));
	}

	generate_draw_line_index();
}

void ElasticPlane::add_spring(const unsigned int index, const double stiffness, const double length, const pair<unsigned int, unsigned int> &between){
	assert(index < springs_.size());

	springs_[index].stiff_ = stiffness;
	springs_[index].length_ = length;
	springs_[index].between_ = between;
}

void ElasticPlane::generate_draw_line(){
	draw_line_.resize(springs_.size()*6, 1);

	for(unsigned int i=0; i<springs_.size(); i++){
		draw_line_.segment<3>(i*6) = points_position_.segment<3>(springs_[i].between_.first*3);
		draw_line_.segment<3>(i*6+3) = points_position_.segment<3>(springs_[i].between_.second*3);
	}
}

void ElasticPlane::generate_draw_line_2D(){
	draw_line_2D_.resize(springs_.size()*6, 1);

	for(unsigned int i=0; i<springs_.size(); i++){
		Vector3d point = Vector3d::Zero();
		point.segment<2>(0) = points_position_2D_.segment<2>(springs_[i].between_.first*2);
		draw_line_2D_.segment<3>(i*6) = point;
		point.segment<2>(0) = points_position_2D_.segment<2>(springs_[i].between_.second*2);
		draw_line_2D_.segment<3>(i*6+3) = point;
	}
}