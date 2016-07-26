#ifndef _FEM_SYSTEM_H_
#define _FEM_SYSTEM_H_

#include <Eigen/Dense>
#include <utility>
#include <set>

#include "../algorithm/point_potential_energy.h"
#include "../algorithm/tetrahedron_potential_energy.h"
#include "abstract_system.h"
#include "point.h"
#include "tetrahedron.h"

namespace zyclincoln{
  class LinearFEMSystem : public AbstractSystem{
  public:
    LinearFEMSystem(const unsigned int points_num, const unsigned int tetrahedrons_num);

    virtual double time_step_ms();
    virtual Eigen::MatrixXd mass_matrix();
    virtual Eigen::VectorXd delta_potential_energy_vector();
    virtual Eigen::MatrixXd delta_delta_potential_energy_matrix();
    virtual Eigen::VectorXd velocity_vector();
    virtual Eigen::VectorXd position_vector();
    virtual void update_velocity_vector(const Eigen::VectorXd &velocity);
    virtual void update_position_vector(const Eigen::VectorXd &position);

    void set_time_step_ms(const double time_step_ms);
    void set_point_potential_energy_calculator(const PointPotentialEnergyCalculator &calculator);
    void set_tetrahedron_potential_energy_calculator(const TetrahedronPotentialEnergyCalculator &calculator);
    Eigen::VectorXd draw_line();
    void add_static_points(std::vector<unsigned int> index_of_static_points);
    void add_points(const std::vector<Point> &points);
    void add_tetrahedrons(const std::vector<Tetrahedron> &tetrahedrons);
    void update_draw_line();
    std::vector<Point> points();
    std::vector<Tetrahedron> tetrahedrons();
  
  private:
    const unsigned int ctetrahedrons_;
    const unsigned int cpoints_;
    double time_step_ms_;
    std::vector<Tetrahedron> tetrahedrons_;
    std::vector<Point> points_;

    Eigen::MatrixXd aux_mass_matrix_;
    Eigen::VectorXd aux_points_velocity_vector_;
    Eigen::VectorXd aux_points_position_vector_;

    Eigen::VectorXd draw_line_;

    std::set<unsigned int> static_points_;

    PointPotentialEnergyCalculator point_potential_energy_calculator_;
    TetrahedronPotentialEnergyCalculator tetrahedron_potential_energy_calculator_;
  };

  inline double LinearFEMSystem::time_step_ms(){
    return time_step_ms_;
  }

  inline Eigen::MatrixXd LinearFEMSystem::mass_matrix(){
    return aux_mass_matrix_;
  }

  inline Eigen::VectorXd LinearFEMSystem::position_vector(){
    return aux_points_position_vector_;
  }

  inline Eigen::VectorXd LinearFEMSystem::velocity_vector(){
    return aux_points_velocity_vector_;
  }

  inline void LinearFEMSystem::set_time_step_ms(const double time_step_ms){
    time_step_ms_ = time_step_ms;
  }

  inline void LinearFEMSystem::set_point_potential_energy_calculator(const PointPotentialEnergyCalculator &calculator){
    point_potential_energy_calculator_ = calculator;
  }

  inline void LinearFEMSystem::set_tetrahedron_potential_energy_calculator(const TetrahedronPotentialEnergyCalculator &calculator){
    tetrahedron_potential_energy_calculator_ = calculator;
  }

  inline Eigen::VectorXd LinearFEMSystem::draw_line(){
    return draw_line_;
  }

  inline std::vector<Point> LinearFEMSystem::points(){
    return points_;
  }

  inline std::vector<Tetrahedron> LinearFEMSystem::tetrahedrons(){
    return tetrahedrons_;
  }

}



#endif