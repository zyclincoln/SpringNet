#ifndef _ELASTIC_SYSTEM_H__
#define _ELASTIC_SYSTEM_H__

#include <vector>
#include <set>
#include <utility>
#include <Eigen/Dense>

#include "../algorithm/spring_potential_energy.h"
#include "../algorithm/point_potential_energy.h"
#include "matrix_op.h"
#include "abstract_system.h"
#include "spring.h"
#include "point.h"

namespace zyclincoln{

	class ElasticSystem : public AbstractSystem{
	public:
		ElasticSystem(const unsigned int points_num, const unsigned int springs_num);
		
		virtual double time_step_ms();
		virtual Eigen::MatrixXd& mass_matrix();
    virtual void delta_potential_energy_vector(Eigen::VectorXd &delta);
    virtual void delta_delta_potential_energy_matrix(Eigen::MatrixXd &delta);
    virtual Eigen::VectorXd& velocity_vector();
    virtual Eigen::VectorXd& position_vector();
    virtual void update_velocity_vector(const Eigen::VectorXd &velocity);
    virtual void update_position_vector(const Eigen::VectorXd &position);

		void set_time_step_ms(const double time_step_ms);
		void set_points_potential_energy_calculator(const PointPotentialEnergyCalculator &calculator);
		void set_springs_potential_energy_calculator(const SpringPotentialEnergyCalculator &calculator);
		Eigen::VectorXd draw_line();
		void add_static_points(std::vector<unsigned int> index_of_static_points);
		void add_points(const std::vector<Point> &points);
		void add_springs(const std::vector<Spring> &springs);
		void update_draw_line();

	private:
		const unsigned int csprings_;
		const unsigned int cpoints_;
		double time_step_ms_;
		std::vector<Spring> springs_;
		std::vector<Point> points_;

		Eigen::MatrixXd aux_mass_matrix_;
		Eigen::VectorXd aux_points_velocity_vector_;
		Eigen::VectorXd aux_points_position_vector_;

		Eigen::VectorXd draw_line_;

		std::set<unsigned int> static_points_;
		std::set<unsigned int> static_lines_;

		PointPotentialEnergyCalculator points_potential_energy_calculator_;
		SpringPotentialEnergyCalculator springs_potential_energy_calculator_;
	};

	inline double ElasticSystem::time_step_ms(){
		return time_step_ms_;
	}

	inline Eigen::MatrixXd& ElasticSystem::mass_matrix(){
		Eigen::MatrixXd shrinked_mass;
		ShrinkMatrix(aux_mass_matrix_, static_lines_, static_lines_, shrinked_mass);

		return shrinked_mass;
	}

	inline Eigen::VectorXd& ElasticSystem::position_vector(){
		Eigen::VectorXd shrinked_position;
		ShrinkColVector(aux_points_position_vector_, static_lines_, shrinked_position);

		return shrinked_position;
	}

	inline Eigen::VectorXd& ElasticSystem::velocity_vector(){
		Eigen::VectorXd shrinked_velocity;
		ShrinkColVector(aux_points_velocity_vector_, static_lines_, shrinked_velocity);

		return shrinked_velocity;
	}

	inline void ElasticSystem::set_time_step_ms(const double time_step_ms){
		time_step_ms_ = time_step_ms;
	}

	inline void ElasticSystem::set_points_potential_energy_calculator(const PointPotentialEnergyCalculator &calculator){
		points_potential_energy_calculator_ = calculator;
	}

	inline void ElasticSystem::set_springs_potential_energy_calculator(const SpringPotentialEnergyCalculator &calculator){
		springs_potential_energy_calculator_ = calculator;
	}

	inline Eigen::VectorXd ElasticSystem::draw_line(){
		return draw_line_;
	}

}

#endif