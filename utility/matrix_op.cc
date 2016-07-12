#include <iostream>
#include "matrix_op.h"

using namespace std;
using namespace Eigen;
using namespace zyclincoln;

using index = unsigned int;

void ShrinkMatrix(const MatrixXd &matrix, const set<index> &remove_row, const set<index> &remove_col, MatrixXd &shrinked_matrix){
	shrinked_matrix = MatrixXd::Zero(matrix.rows() - remove_row.size(), matrix.cols() - remove_col.size());
	
	index original_row = 0, original_col = 0;
	for(index row = 0; row < shrinked_matrix.rows(); row++, original_row++){
		while(remove_row.find(original_row) != remove_row.end()){
			original_row++;
		}
		original_col = 0;
		for(index col = 0; col < shrinked_matrix.cols(); col++, original_col++){
			while(remove_col.find(original_col) != remove_col.end()){
				original_col++;
			}
			shrinked_matrix(row, col) = matrix(original_row, original_col);
		}
	}
}

void ShrinkColvector(const VectorXd &vector, const set<index> &remove_row, VectorXd &shrinked_vector){
	shrinked_vector = VectorXd::Zero(vector.rows() - remove_row.size(), 1);

	index original_row = 0;
	for(index row = 0; row < shrinked_vector.rows(); row++, original_row++){
		while(remove_row.find(original_row)	!= remove_row.end()){
			original_row++;
		}
		shrinked_vector(row, 0) = vector(original_row, 0);
	}
}

void MapToOriginalColVector(const VectorXd &shrinked_vector, const set<index> &remove_row, VectorXd &original_vector){
	original_vector.resize(shrinked_vector.rows() + remove_row.size());
	original_vector.setZero();
	
	index original_row = 0;

	for(index row = 0; row < shrinked_vector.rows(); row++, original_row++){
		while(remove_row.find(original_row) != remove_row.end()){
			original_row++;
		}
		original_vector(original_row, 0) = shrinked_vector(row, 0);
	}
}