#ifndef _EIGEN_UTIL_H_
#define _EIGEN_UTIL_H_

#include <Eigen/Dense>
#include <iostream>
#include <set>

namespace zyclincoln{

  void get_top_k_eigen(Eigen::VectorXd& eigen_value, Eigen::MatrixXd& eigen_vector, 
    Eigen::VectorXd& reduce_value, Eigen::MatrixXd& reduce_vector, unsigned int top_k, unsigned int end_k){

    reduce_value.resize(end_k - top_k, 1);
    reduce_vector.resize(eigen_vector.cols(), end_k - top_k);

    std::set<unsigned int> selected;
    for(unsigned int i = 0; i < end_k; i++){
      double min = 1000000;
      unsigned int min_index = -1;
      for(unsigned int j = 0; j < eigen_value.rows(); j++){
        if(eigen_value(j, 0) < min && selected.find(j) == selected.end()){
          min = eigen_value(j, 0);
          min_index = j;
        }
      }

      if(min_index == -1){
        std::cerr << "[ERROR] cannot find eigen value less than 1000000." << std::endl;
      }
      else{
        selected.insert(min_index);
        if(i >= top_k){
          reduce_value(i-top_k, 0) = min;
          reduce_vector.col(i-top_k) = eigen_vector.col(min_index);
        }
      }

    }

  }
}

#endif