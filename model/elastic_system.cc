#include <assert.h>
#include <math.h>
#include <iostream>
#include "elastic_system.h"
#include "../utility/matrix_op.h"

using namespace Eigen;
using namespace std;
using namespace zyclincoln;

using index = unsigned int;

ElasticSystem::ElasticSystem(const unsigned int points_num, const unsigned int springs_num):
	csprings_(springs_num),
	cpoints_(points_num){
	aux_mass_matrix_ = Eigen::MatrixXd::Zero(cpoints_*3, cpoints_*3);
	aux_points_position_vector_ = Eigen::VectorXd::Zero(cpoints_*3, 1);
	aux_points_velocity_vector_ = Eigen::VectorXd::Zero(cpoints_*3, 1);

	draw_line_ = Eigen::VectorXd::Zero(csprings_*2*3, 1);
}

VectorXd ElasticSystem::delta_potential_energy_vector(){
	VectorXd delta_vector = VectorXd::Zero(cpoints_*3);
	//points potential energy
	for(index i = 0; i < points_.size(); i++){
		delta_vector.segment<3>(3*i) += points_potential_energy_calculator_.calculate_delta_potential_energy(points_[i]);
	}

	//string potential energy
	for(index i = 0; i < strings_.size(); i++){
		VectorXd delta = strings_potential_energy_calculator_.calculate_delta_potential_energy(strings_[i], points_[strings_[i].between.first], points_[strings_[i].between.second]);
		delta_vector.segment<3>(3*strings_[i].between_.first) += delta.segment<3>(0);
		delta_vector.segment<3>(3*strings_[i].between_.second) += delta.segment<3>(3);
	}

	VectorXd shrinked_vector;
	ShrinkColVector(delta_vector, static_lines_, shrinked_vector);
	return shrinked_vector;
}

MatrixXd ElasticSystem::delta_delta_potential_energy_matrix(){
	MatrixXd delta_delta_matrix = MatrixXd::Zero(cpoints_*3, cpoints_*3);
	//points potential energy -- nothing to do

	//string potential energy
	for(index i=0; i < strings_.size(); i++){
		MatrixXd delta = strings_potential_energy_calculator_.calculate_delta_delta_potential_energy(strings_[i], points_[strings_[i].between.first], points_[strings_[i].between.second]);
		delta_delta_matrix.block(3*strings_[i].between_.first, 3*strings_[i].between_.first, 3, 3) = delta.block(0, 0, 3, 3);
		delta_delta_matrix.block(3*strings_[i].between_.first, 3*strings_[i].between_.second, 3, 3) = delta.block(0, 3, 3, 3);
		delta_delta_matrix.block(3*strings_[i].between_.second, 3*strings_[i].between_.first, 3, 3) = delta.block(3, 0, 3, 3);
		delta_delta_matrix.block(3*strings_[i].between_.second, 3*strings_[i].between_.second, 3, 3) = delta.block(3, 3, 3, 3);	
	}

	MatrixXd shrinked_matrix;
	ShrinkMatirx(delta_delta_matrix, static_lines_, static_lines_, shrinked_matrix);

	return shrinked_matrix;
}

void ElasticSystem::update_velocity_vector(const VectorXd &velocity){
	assert(velocity.rows() == points_.size()*3);
	
	VectorXd original_velocity = MapToOriginalColVector(velocity, static_lines_, original_velocity);

	aux_points_velocity_vector_ = original_velocity;
	for(index i = 0; i < points_.size(); i++){
		points_[i].speed_ = original_velocity.segment<3>(i*3);
	}
}

void ElasticSystem::update_position_vector(const VectorXd &position){
	assert(position.rows() == points_.size()*3);

	VectorXd original_position = MapToOriginalColVector(position, static_lines_, original_position);

	aux_points_position_vector_ = original_position;
	for(index i = 0; i < points_.size(); i++){
		points_[i].position_ = original_position.segment<3>(i*3);
	}
}

void ElasticSystem::add_static_points(std::vector<index> index_of_static_points){
	for(index i = 0; i < index_of_static_points.size(); i++){
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

	for(index i = 0; i < points.size(); i++){
		aux_mass_matrix_(current_size + i, current_size + i) = points[i].mass_;
	}
}

void ElasticSystem::add_springs(const std::vector<Spring> &springs){
	assert(springs.size() + springs_.size() <= cstrings_);

	springs_.insert(springs_.end(), springs.begin(), springs.end());
}

void ElasticSystem::update_draw_line(){
	draw_line_ = VectorXd::Zero(springs_.size()*6, 1);

	for(index i=0; i<springs_.size(); i++){
		draw_line_.segment<3>(i*6) = aux_points_position_vector_.segment<3>(springs_[i].between_.first*3);
		draw_line_.segment<3>(i*6+3) = aux_points_position_vector_.segment<3>(springs_[i].between_.second*3);
	}
}