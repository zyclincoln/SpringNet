#ifndef _ELASTIC_PLANE_H__
#define _ELASTIC_PLANE_H__

#include <vector>
#include <set>
#include <utility>
#include <Eigen/Dense>

struct Spring{
public:
	Spring(double stiff, double length, std::pair<unsigned int, unsigned int> between);
	double stiff_;
	double length_;
	std::pair<unsigned int, unsigned int> between_;
};

class ElasticPlane{
public:
	ElasticPlane(const unsigned int length, const unsigned int width, const unsigned int springs);
	// ElasticPlane(const unsigned int length, const unsigned int width);
	ElasticPlane(const unsigned int points_num, const unsigned int springs_num);
	unsigned int clength();
	unsigned int cwidth();
	double time_step_ms();
	Eigen::MatrixXd mass_matrix();
	Eigen::VectorXd points_position();
	Eigen::VectorXd points_speed();
	Eigen::MatrixXd stiffness_matrix();
	Eigen::VectorXd draw_line_index();
	Eigen::VectorXd draw_line();
	void set_time_step_ms(const double time_step_ms);
	void next_frame();
	void add_static_points(std::vector<unsigned int> index_of_static_points);
	void add_points(const std::vector<double> &mass, const std::vector<Eigen::Vector3d> &position, const std::vector<Eigen::Vector3d> &speed);
	void add_point(const unsigned int index, const double mass, const Eigen::Vector3d &position, const Eigen::Vector3d &speed);
	void add_springs(const std::vector<double> &stiffness, const std::vector<double> &length, const std::vector<std::pair<unsigned int, unsigned int> > &between);
	void add_spring(const unsigned int index, const double stiffness, const double length, const std::pair<unsigned int, unsigned int> &between);
	void generate_draw_line();

private:
	void generate_draw_line_index();
	void default_spring_layout();
	void setup_matrix();
	double time_step_ms_;
	const unsigned int clength_;
	const unsigned int cwidth_;
	const unsigned int csprings_;
	Eigen::MatrixXd mass_matrix_;
	Eigen::VectorXd points_position_;
	Eigen::VectorXd points_speed_;
	Eigen::MatrixXd stiffness_matrix_;
	Eigen::VectorXd draw_line_index_;
	Eigen::VectorXd draw_line_;
	std::vector<Spring> springs_;
	std::set<unsigned int> static_points_;
	unsigned int cpoints_;
};

inline Spring::Spring(double stiff, double length, std::pair<unsigned int, unsigned int> between):
	stiff_(stiff),
	length_(length),
	between_(between){

}

inline ElasticPlane::ElasticPlane(const unsigned int length, const unsigned int width, const unsigned springs):
	clength_(length),
	cwidth_(width),
	csprings_(springs),
	cpoints_((length+1)*(width+1)){
	setup_matrix();
}

// inline ElasticPlane::ElasticPlane(const unsigned int length, const unsigned int width):
// 	clength_(length),
// 	cwidth_(width),
// 	csprings_(3*length*width+length+width),
// 	cpoints_((length+1)*(width+1)){
// 	setup_matrix();
// 	default_spring_layout();
// }

inline ElasticPlane::ElasticPlane(const unsigned int points_num, const unsigned int springs_num):
	clength_(0),
	cwidth_(0),
	csprings_(springs_num),
	cpoints_(points_num){
	mass_matrix_.resize(cpoints_*3, cpoints_*3);
	stiffness_matrix_.resize(cpoints_*3, cpoints_*3);
	points_position_.resize(cpoints_*3, 1);
	points_speed_.resize(cpoints_*3, 1);

	mass_matrix_.setZero();
	stiffness_matrix_.setZero();
	points_speed_.setZero();
	points_position_.setZero();
}

inline unsigned int ElasticPlane::clength(){
	return clength_;
}

inline unsigned int ElasticPlane::cwidth(){
	return cwidth_;
}

inline double ElasticPlane::time_step_ms(){
	return time_step_ms_;
}

inline Eigen::MatrixXd ElasticPlane::mass_matrix(){
	return mass_matrix_;
}

inline Eigen::VectorXd ElasticPlane::points_position(){
	return points_position_;
}

inline Eigen::VectorXd ElasticPlane::points_speed(){
	return points_speed_;
}

inline Eigen::MatrixXd ElasticPlane::stiffness_matrix(){
	return stiffness_matrix_;
}

inline Eigen::VectorXd ElasticPlane::draw_line_index(){
	return draw_line_index_;
}

inline void ElasticPlane::set_time_step_ms(const double time_step_ms){
	time_step_ms_ = time_step_ms;
}

inline Eigen::VectorXd ElasticPlane::draw_line(){
	return draw_line_;
}

#endif