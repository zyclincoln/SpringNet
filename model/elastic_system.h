#ifndef _ELASTIC_SYSTEM_H__
#define _ELASTIC_SYSTEM_H__

#include <vector>
#include <set>
#include <utility>
#include <Eigen/Dense>

#include "../algorithm/spring_potential_energy.h"
#include "../algorithm/point_potential_energy.h"
#include "abstract_system.h"
#include "spring.h"
#include "point.h"

namespace zuclincoln{

	class ElasticSystem : public AbstractSystem{
	public:
		ElasticSystem(const unsigned int points_num, const unsigned int springs_num);
		
		virtual double time_step_ms();
		virtual Eigen::MatrixXd mass_matrix();
    virtual Eigen::VectorXd delta_potential_energy_vector();
    virtual Eigen::MatrixXd delta_delta_potential_energy_matrix();
    virtual Eigen::VectorXd velocity_vector();
    virtual Eigen::VectorXd position_vector();
    virtual void update_velocity_vector(const Eigen::VectorXd &velocity);
    virtual void update_position_vector(const Eigen::VectorXd &position);

		void set_time_step_ms(const double time_step_ms);
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

		PointPotentialEnergyCalculator *points_potential_energy_calculator_;
		SpringPotentialEnergyCalculator *springs_potential_energy_calculator_;
	};

	inline double ElasticSystem::time_step_ms(){
		return time_step_ms_;
	}

	inline Eigen::MatrixXd ElasticSystem::mass_matrix(){
		return aux_mass_matrix_;
	}

	inline Eigen::VectorXd ElasticSystem::position_vector(){
		return aux_points_position_vector_;
	}

	inline Eigen::VectorXd ElasticSystem::velocity_vector(){
		return aux_points_velocity_vector_;
	}

	inline void ElasticSystem::set_time_step_ms(const double time_step_ms){
		time_step_ms_ = time_step_ms;
	}

	inline Eigen::VectorXd ElasticSystem::draw_line(){
		return draw_line_;
	}

}

#endif