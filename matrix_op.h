#ifndef _MATRIX_OP_H_
#define _MATRIX_OP_H_

#include <Eigen/Dense>
#include <set>
#include <vector>

void shrink_matrix(const Eigen::MatrixXd &matrix, const std::set<unsigned int> &remove_row, const std::set<unsigned int> &remove_col, Eigen::MatrixXd &shrinked_matrix);
void shrink_colvector(const Eigen::VectorXd &vector, const std::set<unsigned int> &remove_row, Eigen::VectorXd &shrinked_vector);
void map_to_original_colvector(const Eigen::VectorXd &shrinked_vector, const std::set<unsigned int> &remove_row, Eigen::VectorXd &original_vector);


#endif