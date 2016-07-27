#include <assert.h>
#include <math.h>
#include <iostream>
#include "fem_system.h"
#include "matrix_op.h"

using namespace Eigen;
using namespace std;
using namespace zyclincoln;

LinearFEMSystem::LinearFEMSystem(const unsigned int points_num, const unsigned int tetrahedrons_num):
  ctetrahedrons_(tetrahedrons_num),
  cpoints_(points_num){
  aux_mass_matrix_ = Eigen::MatrixXd::Zero(cpoints_*3, cpoints_*3);
  aux_points_position_vector_ = Eigen::VectorXd::Zero(cpoints_*3, 1);
  aux_points_velocity_vector_ = Eigen::VectorXd::Zero(cpoints_*3, 1);

  draw_line_ = Eigen::VectorXd::Zero(ctetrahedrons_*6*2*3, 1);
}

VectorXd LinearFEMSystem::delta_potential_energy_vector(){
  VectorXd delta_vector = VectorXd::Zero(cpoints_*3);
  //points potential energy
  for(unsigned int i = 0; i < points_.size(); i++){
    delta_vector.segment<3>(3*i) += point_potential_energy_calculator_.calculate_delta_potential_energy(points_[i]);
  }

  //tetrahedron potential energy
  for(unsigned int i = 0; i < tetrahedrons_.size(); i++){
    tetrahedron_potential_energy_calculator_.PreCompute(tetrahedrons_[i], points_[tetrahedrons_[i].points_index_[0]],
      points_[tetrahedrons_[i].points_index_[1]], points_[tetrahedrons_[i].points_index_[2]], points_[tetrahedrons_[i].points_index_[3]]);
    VectorXd delta = tetrahedron_potential_energy_calculator_.calculate_delta_potential_energy(tetrahedrons_[i]);
    for(unsigned int j = 0; j < 4; j++){
      delta_vector.segment<3>(3*tetrahedrons_[i].points_index_[j]) += delta.segment<3>(j*3);
    }
  }

  return delta_vector;
}

MatrixXd LinearFEMSystem::delta_delta_potential_energy_matrix(){
  if(delta_2_potential_matrix_computed_ == false){
    delta_2_potential_matrix_ = MatrixXd::Zero(cpoints_*3, cpoints_*3);
    //points potential energy -- nothing to do

    //tetrahedron potential energy
    for(unsigned int i = 0; i < tetrahedrons_.size(); i++){
      MatrixXd delta = tetrahedron_potential_energy_calculator_.calculate_delta_delta_potential_energy(tetrahedrons_[i]);
      for(unsigned int j = 0; j < 4; j++){
        for(unsigned int k = 0; k < 4; k++){
          delta_2_potential_matrix_.block(3*tetrahedrons_[i].points_index_[j],3*tetrahedrons_[i].points_index_[k], 3, 3) += delta.block(j*3, k*3, 3, 3);
        }
      }
    }
    delta_2_potential_matrix_computed_ = true;
  }

  return delta_2_potential_matrix_;
}

void LinearFEMSystem::update_velocity_vector(const VectorXd &velocity){
  assert(velocity.rows() == points_.size()*3);

  for(unsigned int i = 0; i < points_.size(); i++){
    if(static_points_.find(i) == static_points_.end()){
      points_[i].speed_ = velocity.segment<3>(i*3);
      aux_points_velocity_vector_.segment<3>(i*3) = velocity.segment<3>(i*3);
    }
  }
}

void LinearFEMSystem::update_position_vector(const VectorXd &position){
  assert(position.rows() == points_.size()*3);

  for(unsigned int i = 0; i < points_.size(); i++){
    if(static_points_.find(i) == static_points_.end()){
      points_[i].position_ = position.segment<3>(i*3);
      aux_points_position_vector_.segment<3>(i*3) = position.segment<3>(i*3);
    }
  }
}

void LinearFEMSystem::add_static_points(std::vector<unsigned int> index_of_static_points){
  for(unsigned int i = 0; i < index_of_static_points.size(); i++){
    static_points_.insert(index_of_static_points[i]);
  }
}

void LinearFEMSystem::add_points(const std::vector<Point> &points){
  assert(points.size() + points_.size() <= cpoints_);

  unsigned int current_size = points_.size();

  points_.insert(points_.end(), points.begin(), points.end());

  for(unsigned int i = 0; i < points.size(); i++){
    aux_mass_matrix_.block((current_size + i)*3, (current_size + i)*3, 3, 3) = MatrixXd::Identity(3, 3)*points[i].mass_;
    aux_points_position_vector_.segment<3>((current_size + i)*3) = points[i].position_;
    aux_points_velocity_vector_.segment<3>((current_size + i)*3) = points[i].speed_;
  }

  cout<<"add "<<points.size()<<" points"<<endl;
}

void LinearFEMSystem::add_tetrahedrons(const std::vector<Tetrahedron> &tetrahedrons){
  assert(tetrahedrons.size() + tetrahedrons_.size() <= ctetrahedrons_);

  unsigned int current_size = tetrahedrons_.size();

  tetrahedrons_.insert(tetrahedrons_.end(), tetrahedrons.begin(), tetrahedrons.end());

  for(unsigned int i = 0; i < tetrahedrons.size(); i++){
    vector<Point> point_set;
    point_set.push_back(points_[tetrahedrons[i].points_index_[0]]);
    point_set.push_back(points_[tetrahedrons[i].points_index_[1]]);
    point_set.push_back(points_[tetrahedrons[i].points_index_[2]]);
    point_set.push_back(points_[tetrahedrons[i].points_index_[3]]);
    tetrahedrons_[current_size + i].PreCompute(point_set);
  }
}

void LinearFEMSystem::update_draw_line(){
  draw_line_ = VectorXd::Zero(tetrahedrons_.size()*6*2*3, 1);
  for(unsigned int i = 0; i<tetrahedrons_.size(); i++){
    draw_line_.segment<3>(i*36) = aux_points_position_vector_.segment<3>(tetrahedrons_[i].points_index_[0]);
    draw_line_.segment<3>(i*36+3) = aux_points_position_vector_.segment<3>(tetrahedrons_[i].points_index_[1]);
    draw_line_.segment<3>(i*36+6) = aux_points_position_vector_.segment<3>(tetrahedrons_[i].points_index_[0]);
    draw_line_.segment<3>(i*36+9) = aux_points_position_vector_.segment<3>(tetrahedrons_[i].points_index_[2]);
    draw_line_.segment<3>(i*36+12) = aux_points_position_vector_.segment<3>(tetrahedrons_[i].points_index_[0]);
    draw_line_.segment<3>(i*36+15) = aux_points_position_vector_.segment<3>(tetrahedrons_[i].points_index_[3]);
    draw_line_.segment<3>(i*36+18) = aux_points_position_vector_.segment<3>(tetrahedrons_[i].points_index_[1]);
    draw_line_.segment<3>(i*36+21) = aux_points_position_vector_.segment<3>(tetrahedrons_[i].points_index_[2]);
    draw_line_.segment<3>(i*36+24) = aux_points_position_vector_.segment<3>(tetrahedrons_[i].points_index_[1]);
    draw_line_.segment<3>(i*36+27) = aux_points_position_vector_.segment<3>(tetrahedrons_[i].points_index_[3]);
    draw_line_.segment<3>(i*36+30) = aux_points_position_vector_.segment<3>(tetrahedrons_[i].points_index_[2]);
    draw_line_.segment<3>(i*36+33) = aux_points_position_vector_.segment<3>(tetrahedrons_[i].points_index_[3]);
  }
}