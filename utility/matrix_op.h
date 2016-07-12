#ifndef _MATRIX_OP_H_
#define _MATRIX_OP_H_

#include <Eigen/Dense>
#include <set>
#include <vector>

namespace zyclincoln{

  void ShrinkMatrix(const Eigen::MatrixXd &matrix, const std::set<unsigned int> &remove_row, const std::set<unsigned int> &remove_col, Eigen::MatrixXd &shrinked_matrix);
  void ShrinkColVector(const Eigen::VectorXd &vector, const std::set<unsigned int> &remove_row, Eigen::VectorXd &shrinked_vector);
  void MapToOriginalColVector(const Eigen::VectorXd &shrinked_vector, const std::set<unsigned int> &remove_row, Eigen::VectorXd &original_vector);

}

#endif