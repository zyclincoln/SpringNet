#include <assert.h>
#include <math.h>
#include <iostream>
#include <Eigen/Sparse>
#include "elastic_plane.h"
#include "matrix_op.h"

using namespace Eigen;
using namespace std;

void ElasticPlane::default_spring_layout(){
	for(unsigned int i = 0; i <= cwidth_; i++){
		for(unsigned int j = 0; j < clength_; j++){
			springs_.push_back(Spring(1, 1, pair<unsigned int, unsigned int>(i*(clength_+1) + j, i*(clength_+1) + j + 1)));
			cout<<i*(clength_+1) + j<<" "<<i*(clength_+1) + j + 1<<endl;
		}
	}

	for(unsigned int i = 0; i < cwidth_; i++){
		for(unsigned int j = 0; j <= clength_; j++){
			springs_.push_back(Spring(1, 1, pair<unsigned int, unsigned int>(i*(clength_+1) + j, (i+1)*(clength_+1) + j)));
			cout<<i*(clength_+1) + j<<" "<<(i+1)*(clength_+1) + j<<endl;
		}
	}

	for(unsigned int i=0; i < cwidth_ ;i++){
		for(unsigned int j = 0; j < clength_; j++){
			springs_.push_back(Spring(1, 1, pair<unsigned int, unsigned int>(i*(clength_+1) + j, (i+1)*(clength_+1) + j + 1)));
			cout<<i*(clength_+1) + j<<" "<<(i+1)*(clength_+1) + j + 1<<endl;
		}
	}

	generate_draw_line_index();
}

void ElasticPlane::generate_draw_line_index(){
	draw_line_index_ = VectorXd::Zero(springs_.size()*2, 1);

	for(unsigned int i = 0; i < springs_.size(); i++){
		draw_line_index_(i*2 + 0, 0) = springs_[i].between_.first;
		draw_line_index_(i*2 + 1, 0) = springs_[i].between_.second;
	}
}

void ElasticPlane::setup_matrix(){
	mass_matrix_.resize((clength_+1)*(cwidth_+1)*3, (clength_+1)*(cwidth_+1)*3);
	stiffness_matrix_.resize((clength_+1)*(cwidth_+1)*3, (clength_+1)*(cwidth_+1)*3);
	points_position_.resize((clength_+1)*(cwidth_+1)*3, 1);
	points_speed_.resize((clength_+1)*(cwidth_+1)*3, 1);

	stiffness_matrix_.setZero();
	mass_matrix_.setZero();
	points_speed_.setZero();
	points_position_.setZero();
}

void ElasticPlane::next_frame(){

  VectorXd xstar = points_position_;

  set<unsigned int> index_need_remove;
  for(unsigned int i=0; i < cpoints_; i++){
    if(static_points_.find(i)!=static_points_.end()){
      index_need_remove.insert(i*3);
      index_need_remove.insert(i*3+1);
      index_need_remove.insert(i*3+2);
    }
  }
  
  for (int iter = 0; iter < 10; ++iter) {
    cout << "iter " << iter << endl;
        stiffness_matrix_.setZero();
	VectorXd force = VectorXd::Zero(cpoints_*3, 1);
	for(unsigned int i = 0; i < springs_.size(); i++){
		unsigned index0 = springs_[i].between_.first;
		unsigned index1 = springs_[i].between_.second;
		Vector3d p0 = xstar.segment<3>(3*index0); //xn.block(index0*3, 0, 3, 1);
		Vector3d p1 = xstar.segment<3>(3*index1); //xn.block(index1*3, 0, 3, 1);

		Matrix3d sub_stiff = Matrix3d::Identity()*(1 - 1/springs_[i].length_);
		sub_stiff += (p1 - p0)*(p1 - p0).transpose()/pow(springs_[i].length_,3);
		// cout<<"===sub_stiff===\n"<<sub_stiff<<endl;
		stiffness_matrix_.block(index0*3, index0*3, 3, 3) += sub_stiff;
		stiffness_matrix_.block(index1*3, index1*3, 3, 3) += sub_stiff;
		stiffness_matrix_.block(index0*3, index1*3, 3, 3) -= sub_stiff;
		stiffness_matrix_.block(index1*3, index0*3, 3, 3) -= sub_stiff;
		// cout<<"===stiffness_matrix===\n"<<stiffness_matrix_<<endl;
		force.block(index0*3, 0, 3, 1) += ((p1 - p0).norm() - springs_[i].length_)*springs_[i].stiff_*(p1 - p0).normalized();
		force.block(index1*3, 0, 3, 1) += ((p1 - p0).norm() - springs_[i].length_)*springs_[i].stiff_*(p0 - p1).normalized();
		// cout<<"===p_force===\n"<<((p1 - p0).norm() - springs_[i].length_)*springs_[i].stiff_*(p1 - p0).normalized()<<endl;
		// cout<<"===p_force===\n"<<((p1 - p0).norm() - springs_[i].length_)*springs_[i].stiff_*(p0 - p1).normalized()<<endl;
		cout<<"===part force===\n"<<index0<<" : "<<((p1 - p0).norm() - springs_[i].length_)*springs_[i].stiff_*(p1 - p0).normalized().transpose()<<endl;
		cout<<index1<<" : "<<((p1 - p0).norm() - springs_[i].length_)*springs_[i].stiff_*(p0 - p1).normalized().transpose()<<endl;
	}

	for(unsigned int i=0; i < cpoints_; i++){
		force.block(i*3, 0, 3, 1) += Vector3d(0, -9.8*0.1, 0);
	}

	cout<<"===total force===\n"<<force<<endl;

	MatrixXd A = MatrixXd::Zero(cpoints_*3, cpoints_*3);
	A = mass_matrix_ + pow(time_step_ms_, 2)*stiffness_matrix_;
	VectorXd b = VectorXd::Zero(cpoints_*3, 1);
        //b = time_step_ms_*(force - stiffness_matrix_*time_step_ms_*points_speed_);
        b = -mass_matrix_*(xstar-points_position_-time_step_ms_*points_speed_)+time_step_ms_*time_step_ms_*force;

	MatrixXd shrinked_A;
	VectorXd shrinked_b;

	shrink_matrix(A, index_need_remove, index_need_remove, shrinked_A);
	shrink_colvector(b, index_need_remove, shrinked_b);
	// cout<<"===force===\n"<<force<<endl;
	// cout<<"====A===\n"<<A<<endl;
	// cout<<"===b===\n"<<b<<endl;

        // Eigen::SparseMatrix<double> AA = A.sparseView();
        // SimplicialCholesky<SparseMatrix<double>> solver;
        // solver.compute(AA);
	VectorXd shrinked_delta_x = shrinked_A.colPivHouseholderQr().solve(shrinked_b);// solver.solve(shrinked_b);
	VectorXd delta_x = VectorXd::Zero(cpoints_*3, 1);
        map_to_original_colvector(shrinked_delta_x, index_need_remove, delta_x);
        xstar += delta_x;
       
        //xn = points_position_+vn*time_step_ms_;
	// cout<<"===delta v===\n"<<delta_v<<endl;

     }

  points_speed_ = (xstar-points_position_)/time_step_ms_;
  points_position_ = xstar; //points_speed_*time_step_ms_;
  // points_speed_ += delta_v;
  //points_position_ += points_speed_*time_step_ms_;
	cout<<"===point position===\n"<<points_position_<<endl;
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

void ElasticPlane::add_point(const unsigned int index, const double mass, const Vector3d &position, const Vector3d &speed){
	assert(index < (clength_+1)*(cwidth_+1));
	
	mass_matrix_.block(index*3, index*3, 3, 3) = Matrix3d::Identity()*mass;
	points_position_.block(index*3, 0, 3, 1) = position;
	points_speed_.block(index*3, 0, 3, 1) = speed;
}

void ElasticPlane::add_springs(const vector<double> &stiffness, const vector<double> &length, const vector<pair<unsigned int, unsigned int> > &between){
	assert(stiffness.size() == length.size());
	assert(stiffness.size() == between.size());

	for(unsigned int i = 0; i < stiffness.size(); i++){
          springs_.push_back(Spring(1e2, /* stiffness[i],*/ length[i], between[i]));
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
		draw_line_.block(i*6, 0, 3, 1) = points_position_.block(springs_[i].between_.first*3, 0, 3, 1);
		draw_line_.block(i*6+3, 0, 3, 1) = points_position_.block(springs_[i].between_.second*3, 0, 3, 1);
	}
}
