#include <sys/time.h>

#include <assert.h>
#include <math.h>
#include <iostream>
#include <utility>
#include <omp.h>
#include "fem_system.h"
#include "../algorithm/spring_potential_energy.h"
#include "spring.h"
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

void LinearFEMSystem::delta_potential_energy_vector(VectorXd &delta_vector){
  delta_vector = VectorXd::Zero(cpoints_*3);
  //points potential energy
  for(unsigned int i = 0; i < points_.size(); i++){
    delta_vector.segment<3>(3*i) += point_potential_energy_calculator_.calculate_delta_potential_energy(points_[i]);
  }

  //tetrahedron potential energy
  for(unsigned int i = 0; i < tetrahedrons_.size(); i++){
    tetrahedron_potential_energy_calculator_.PreCompute(tetrahedrons_[i], points_[tetrahedrons_[i].points_index_[0]],
      points_[tetrahedrons_[i].points_index_[1]], points_[tetrahedrons_[i].points_index_[2]], points_[tetrahedrons_[i].points_index_[3]], corotate_);

    // cout << " point :\n" << points_[tetrahedrons_[i].points_index_[0]].position_.transpose() << endl << points_[tetrahedrons_[i].points_index_[1]].position_.transpose() << endl
    //   << points_[tetrahedrons_[i].points_index_[2]].position_.transpose() << endl << points_[tetrahedrons_[i].points_index_[3]].position_.transpose() << endl;
    // cout << "inverse: \n" << tetrahedrons_[i].p_inverse_ << endl;
    // Matrix<double, 12, 1> delta;
    VectorXd delta;
    delta.resize(12, 1);
    tetrahedron_potential_energy_calculator_.calculate_delta_potential_energy(tetrahedrons_[i], delta, corotate_);
    for(unsigned int j = 0; j < 4; j++){
      delta_vector.segment<3>(3*tetrahedrons_[i].points_index_[j]) += delta.segment<3>(j*3);
    }
  }
    //constraint as string
  SimpleSpringPotentialEnergyCalculator constraint_calculator;
  for(unsigned int i = 0; i < static_points_vector_.size(); i++){
    Spring spring(5000, 0, pair<unsigned int, unsigned int>(0, 1));
    // Vector3d delta;
    // delta.setZero();
    VectorXd delta;
    delta.resize(3, 1);
    constraint_calculator.calculate_delta_potential_energy(spring, points_[static_points_vector_[i]], static_points_original_[i], delta);
    delta_vector.segment<3>(3*static_points_vector_[i]) += delta.segment<3>(0); 
    // cout << "spring " << i << endl << delta.transpose() << endl;
  }
}

void LinearFEMSystem::delta_delta_potential_energy_matrix(MatrixXd &delta_2_potential_matrix){
  timeval time0, time1;
  if(delta_2_potential_matrix_computed_ == false || corotate_ == true){
    delta_2_potential_matrix_ = MatrixXd::Zero(cpoints_*3, cpoints_*3);
    //points potential energy -- nothing to do

    //tetrahedron potential energy

    gettimeofday(&time0, 0);

    for(unsigned int i = 0; i < tetrahedrons_.size(); i++){
      // Matrix<double, 12, 12> delta;
      MatrixXd delta;
      delta.resize(12, 12);
      tetrahedron_potential_energy_calculator_.calculate_delta_delta_potential_energy(tetrahedrons_[i], delta, corotate_);
      for(unsigned int k = 0; k < 4; k++){
        for(unsigned int j = 0; j < 4; j++){
          delta_2_potential_matrix_.block(3*tetrahedrons_[i].points_index_[j],3*tetrahedrons_[i].points_index_[k], 3, 3) += delta.block(j*3, k*3, 3, 3);
        }
      }
    }
    gettimeofday(&time1, 0);

    // cout << "== time to cal del2: " << (1000000*(time1.tv_sec - time0.tv_sec) + time1.tv_usec - time0.tv_usec)/1000 << " ms" << endl;
    delta_2_potential_matrix_computed_ = true;

    //constraint as string
    SimpleSpringPotentialEnergyCalculator constraint_calculator;
    for (unsigned int i = 0; i < static_points_vector_.size(); ++i){
      Spring spring(5000, 0, pair<unsigned int, unsigned int>(0, 1));
      // Matrix3d delta;
      // delta.setZero();
      MatrixXd delta;
      delta.resize(3, 3);
      constraint_calculator.calculate_delta_delta_potential_energy(spring, points_[static_points_vector_[i]], static_points_original_[i], delta);
      delta_2_potential_matrix_.block(3*static_points_vector_[i], 3*static_points_vector_[i], 3, 3) += delta.block(0, 0, 3, 3); 
    }
  }
  delta_2_potential_matrix = delta_2_potential_matrix_;
}

void LinearFEMSystem::update_velocity_vector(const VectorXd &velocity){
  assert(velocity.rows() == points_.size()*3);

  for(unsigned int i = 0; i < points_.size(); i++){
    // if(static_points_.find(i) == static_points_.end()){
      points_[i].speed_ = velocity.segment<3>(i*3);
      aux_points_velocity_vector_.segment<3>(i*3) = velocity.segment<3>(i*3);
    // }
  }
}

void LinearFEMSystem::update_position_vector(const VectorXd &position){
  assert(position.rows() == points_.size()*3);

  for(unsigned int i = 0; i < points_.size(); i++){
    // if(static_points_.find(i) == static_points_.end()){
      points_[i].position_ = position.segment<3>(i*3);
      aux_points_position_vector_.segment<3>(i*3) = position.segment<3>(i*3);
    // }
  }
}

void LinearFEMSystem::add_static_points(std::vector<unsigned int> index_of_static_points){
  for(unsigned int i = 0; i < index_of_static_points.size(); i++){
    static_points_.insert(index_of_static_points[i]);
    static_points_vector_.push_back(index_of_static_points[i]);
    static_points_original_.push_back(points_[index_of_static_points[i]]);
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

  // cout<<"add "<<points.size()<<" points"<<endl;
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
    // cout << "index: " << i << endl;
    // cout << "point: " << points_[tetrahedrons[i].points_index_[0]].position_.transpose() << endl;
    // cout << "point: " << points_[tetrahedrons[i].points_index_[1]].position_.transpose() << endl;
    // cout << "point: " << points_[tetrahedrons[i].points_index_[2]].position_.transpose() << endl;
    // cout << "point: " << points_[tetrahedrons[i].points_index_[3]].position_.transpose() << endl;
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

double LinearFEMSystem::total_energy(){
  double point_potential_energy = 0;
  double point_kinetic_energy = 0;
  for(unsigned int i = 0; i < points_.size(); i++){
    point_potential_energy += point_potential_energy_calculator_.calculate_potential_energy(points_[i]);
    point_kinetic_energy += 0.5*points_[i].mass_*pow(points_[i].speed_.norm(), 2);
  }


  double string_potential_energy = 0;
  for(unsigned int i = 0; i < static_points_vector_.size(); i++){
    SpringPotentialEnergyCalculator calculator;
    Spring spring(100, 0, pair<unsigned int, unsigned int>(0, 1));
    string_potential_energy += calculator.calculate_potential_energy(spring, points_[static_points_vector_[i]], static_points_original_[i]);
  }

  double tetrahedron_potential_energy = 0;
  for(unsigned int i = 0; i < tetrahedrons_.size(); i++){
    tetrahedron_potential_energy += tetrahedron_potential_energy_calculator_.calculate_potential_energy(tetrahedrons_[i]);
  }

  // cout << "=== energy report begin===" << endl;
  // cout << "point potential energy : " << point_potential_energy << endl;
  // cout << "point kinetic energy : " << point_kinetic_energy << endl;
  // cout << "string potential energy : " << string_potential_energy << endl; 
  // cout << "tetrahedron potential energy : " << tetrahedron_potential_energy << endl;
  // cout << "total energy : " << point_kinetic_energy + point_potential_energy+ string_potential_energy + tetrahedron_potential_energy << endl;
  // cout << "=== energy report end ===" << endl;
  return point_kinetic_energy + point_potential_energy+ string_potential_energy + tetrahedron_potential_energy;
}