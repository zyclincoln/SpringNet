#ifndef _MATRIX_OP_H_
#define _MATRIX_OP_H_

#include <Eigen/Dense>
#include <set>
#include <vector>

Eigen::MatrixXd shrink_matrix(const Eigen::MatrixXd &matrix, const std::set<unsigned int> &update_row, const std::set<unsigned int> &update_col);
Eigen::VectorXd shrink_colvector(const Eigen::VectorXd &vector, const std::set<unsigned int> &update_row);
void map_to_original_colvector(const Eigen::VectorXd &shrinked_vector, const std::set<unsigned int> &update_row, Eigen::VectorXd &original_vector);


#endif