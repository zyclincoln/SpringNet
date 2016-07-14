#include <Eigen/Dense>
#include "spring.h"

using namespace std;
using namespace zyclincoln;

Spring::Spring(const double stiff, const double length, const pair<unsigned int, unsigned int> between):
  stiff_(stiff),
  length_(length),
  between_(between){

}