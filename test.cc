#include <Eigen/Dense>
#include <iostream>
#include <Eigen/SVD>

using namespace std;
using namespace Eigen;

int main(){
  Matrix3d temp;
  temp = Matrix3d::Random();

  JacobiSVD<MatrixXd> svd(temp, ComputeThinU | ComputeThinV);
  cout << "===singular value===" << endl << svd.singularValues() << endl;
  cout << "===matrix U===" << endl << svd.matrixU() << endl;
  cout << "===matrix V===" << endl << svd.matrixV() << endl;

  Matrix3d singular;
  singular = Matrix3d::Zero();
  singular(0, 0) = svd.singularValues()(0, 0);
  singular(1, 1) = svd.singularValues()(1, 0);
  singular(2, 2) = svd.singularValues()(2, 0);
  
  cout << "===result===" << endl << svd.matrixU()*singular*svd.matrixV().transpose() - temp << endl;
}