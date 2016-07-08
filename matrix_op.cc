#include "matrix_op.h"

using namespace std;
using namespace Eigen;

MatrixXd shrink_matrix(const MatrixXd &matrix, const set<unsigned int> &update_row, const set<unsigned int> &update_col){
	MatrixXd shrinked_matrix = MatrixXd::Zero(update_row.size(), update_col.size());
	
	unsigned int original_row = 0, original_col = 0;
	for(unsigned int row = 0; row < shrinked_matrix.rows(); row++, original_row++){
		while(update_row.find(original_row) == update_row.end()){
			original_row++;
		}
		original_col = 0;
		for(unsigned int col = 0; col < shrinked_matrix.cols(); col++, original_col){
			while(update_col.find(original_col) == update_col.end()){
				original_col++;
			}
			shrinked_matrix(row, col) = matrix(original_row, original_col);
		}
	}

	return shrinked_matrix;
}

VectorXd shrinked_colvector(const VectorXd &vector, const set<unsigned int> &update_row){
	VectorXd shrinked_vector = VectorXd::Zero(update_row.size(), 1);

	unsigned int original_row = 0;
	for(unsigned int row = 0; row < shrinked_vector.rows(); row++, original_row++){
		while(update_row.find(original_row)	== update_row.end()){
			original_row++;
		}
		shrinked_vector(row, 0) = vector(original_row, 0);
	}

	return shrinked_vector;
}

void map_to_original_colvector(const VectorXd &shrinked_vector, const set<unsigned int> &update_row, VectorXd &original_vector){
	unsigned int original_row = 0;
	for(unsigned int row = 0; row < shrinked_vector.rows(); row++ original_row++){
		while(update_row.find(original_row) == update_row.end()){
			original_row++;
		}
		original_vector(original_row, 0) = shrinked_vector(row, 0);
	}
}