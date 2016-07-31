#include <assert.h>
#include <math.h>
#include <iostream>
#include "elastic_system.h"
#include "matrix_op.h"

using namespace Eigen;
using namespace std;
using namespace zyclincoln;

ElasticSystem::ElasticSystem(const unsigned int points_num, const unsigned int springs_num):
	csprings_(springs_num),
	cpoints_(points_num){
	aux_mass_matrix_ = Eigen::MatrixXd::Zero(cpoints_*3, cpoints_*3);
	aux_points_position_vector_ = Eigen::VectorXd::Zero(cpoints_*3, 1);
	aux_points_velocity_vector_ = Eigen::VectorXd::Zero(cpoints_*3, 1);

	draw_line_ = Eigen::VectorXd::Zero(csprings_*2*3, 1);
}

void ElasticSystem::delta_potential_energy_vector(VectorXd &delta){
	VectorXd delta_vector = VectorXd::Zero(cpoints_*3);
	//points potential energy
	for(unsigned int i = 0; i < points_.size(); i++){
		delta_vector.segment<3>(3*i) += points_potential_energy_calculator_.calculate_delta_potential_energy(points_[i]);
	}

	//string potential energy
	for(unsigned int i = 0; i < springs_.size(); i++){
		VectorXd sub_delta;
		springs_potential_energy_calculator_.calculate_delta_potential_energy(springs_[i], points_[springs_[i].between_.first], points_[springs_[i].between_.second], sub_delta);
		delta_vector.segment<3>(3*springs_[i].between_.first) += sub_delta.segment<3>(0);
		delta_vector.segment<3>(3*springs_[i].between_.second) += sub_delta.segment<3>(3);
	}

	VectorXd shrinked_vector;
	ShrinkColVector(delta_vector, static_lines_, shrinked_vector);
	delta = shrinked_vector;
}

void ElasticSystem::delta_delta_potential_energy_matrix(MatrixXd &delta){
	MatrixXd delta_delta_matrix = MatrixXd::Zero(cpoints_*3, cpoints_*3);
	//points potential energy -- nothing to do

	//string potential energy
	for(unsigned int i=0; i < springs_.size(); i++){
		// Matrix<double, 6, 6> delta;
		MatrixXd delta;
		delta.resize(6, 6);
		springs_potential_energy_calculator_.calculate_delta_delta_potential_energy(springs_[i], points_[springs_[i].between_.first], points_[springs_[i].between_.second], delta);
		delta_delta_matrix.block(3*springs_[i].between_.first, 3*springs_[i].between_.first, 3, 3) += delta.block(0, 0, 3, 3);
		delta_delta_matrix.block(3*springs_[i].between_.first, 3*springs_[i].between_.second, 3, 3) += delta.block(0, 3, 3, 3);
		delta_delta_matrix.block(3*springs_[i].between_.second, 3*springs_[i].between_.first, 3, 3) += delta.block(3, 0, 3, 3);
		delta_delta_matrix.block(3*springs_[i].between_.second, 3*springs_[i].between_.second, 3, 3) += delta.block(3, 3, 3, 3);	
	}

	MatrixXd shrinked_matrix;
	ShrinkMatrix(delta_delta_matrix, static_lines_, static_lines_, shrinked_matrix);

	return shrinked_matrix;
}

void ElasticSystem::update_velocity_vector(const VectorXd &velocity){
	VectorXd original_velocity;
	MapToOriginalColVector(velocity, static_lines_, original_velocity);
	assert(original_velocity.rows() == points_.size()*3);

	for(unsigned int i = 0; i < points_.size(); i++){
		if(static_points_.find(i) == static_points_.end()){
			points_[i].speed_ = original_velocity.segment<3>(i*3);
			aux_points_velocity_vector_.segment<3>(i*3) = original_velocity.segment<3>(i*3);
		}
	}
}

void ElasticSystem::update_position_vector(const VectorXd &position){
	VectorXd original_position;
	MapToOriginalColVector(position, static_lines_, original_position);
	assert(original_position.rows() == points_.size()*3);

	for(unsigned int i = 0; i < points_.size(); i++){
		if(static_points_.find(i) == static_points_.end()){
			points_[i].position_ = original_position.segment<3>(i*3);
			aux_points_position_vector_.segment<3>(i*3) = original_position.segment<3>(i*3);
		}
	}
}

void ElasticSystem::add_static_points(std::vector<unsigned int> index_of_static_points){
	for(unsigned int i = 0; i < index_of_static_points.size(); i++){
		static_points_.insert(index_of_static_points[i]);
		static_lines_.insert(index_of_static_points[i]*3);
		static_lines_.insert(index_of_static_points[i]*3+1);
		static_lines_.insert(index_of_static_points[i]*3+2);
	}
}

void ElasticSystem::add_points(const std::vector<Point> &points){
	assert(points.size() + points_.size() <= cpoints_);

	unsigned int current_size = points_.size();

	points_.insert(points_.end(), points.begin(), points.end());

	for(unsigned int i = 0; i < points.size(); i++){
		aux_mass_matrix_.block((current_size + i)*3, (current_size + i)*3, 3, 3) = MatrixXd::Identity(3, 3)*points[i].mass_;
		aux_points_position_vector_.segment<3>((current_size + i)*3) = points[i].position_;
		aux_points_velocity_vector_.segment<3>((current_size + i)*3) = points[i].speed_;
	}
}

void ElasticSystem::add_springs(const std::vector<Spring> &springs){
	assert(springs.size() + springs_.size() <= csprings_);

	springs_.insert(springs_.end(), springs.begin(), springs.end());
}

void ElasticSystem::update_draw_line(){
	draw_line_ = VectorXd::Zero(springs_.size()*6, 1);

	for(unsigned int i=0; i<springs_.size(); i++){
		draw_line_.segment<3>(i*6) = aux_points_position_vector_.segment<3>(springs_[i].between_.first*3);
		draw_line_.segment<3>(i*6+3) = aux_points_position_vector_.segment<3>(springs_[i].between_.second*3);
	}
}