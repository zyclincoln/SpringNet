#include <iostream>
#include "matrix_op.h"

using namespace std;
using namespace Eigen;

void shrink_matrix(const MatrixXd &matrix, const set<unsigned int> &remove_row, const set<unsigned int> &remove_col, MatrixXd &shrinked_matrix){
	shrinked_matrix = MatrixXd::Zero(matrix.rows() - remove_row.size(), matrix.cols() - remove_col.size());
	
	unsigned int original_row = 0, original_col = 0;
	for(unsigned int row = 0; row < shrinked_matrix.rows(); row++, original_row++){
		while(remove_row.find(original_row) != remove_row.end()){
			original_row++;
		}
		original_col = 0;
		for(unsigned int col = 0; col < shrinked_matrix.cols(); col++, original_col++){
			while(remove_col.find(original_col) != remove_col.end()){
				original_col++;
			}
			shrinked_matrix(row, col) = matrix(original_row, original_col);
		}
	}
}

void shrink_colvector(const VectorXd &vector, const set<unsigned int> &remove_row, VectorXd &shrinked_vector){
	shrinked_vector = VectorXd::Zero(vector.rows() - remove_row.size(), 1);

	unsigned int original_row = 0;
	for(unsigned int row = 0; row < shrinked_vector.rows(); row++, original_row++){
		while(remove_row.find(original_row)	!= remove_row.end()){
			original_row++;
		}
		shrinked_vector(row, 0) = vector(original_row, 0);
	}
}

void map_to_original_colvector(const VectorXd &shrinked_vector, const set<unsigned int> &remove_row, VectorXd &original_vector){
	unsigned int original_row = 0;
	// cout<<"====shrinked_vector====\n"<<shrinked_vector.transpose()<<endl;
	for(unsigned int row = 0; row < shrinked_vector.rows(); row++, original_row++){
		while(remove_row.find(original_row) != remove_row.end()){
			original_row++;
		}
		original_vector(original_row, 0) = shrinked_vector(row, 0);
		// cout<<"====map process===\n"<<original_vector.transpose()<<endl;
	}
}