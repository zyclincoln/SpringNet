#include <Eigen/Dense>
#include "spring.h"

using namespace zyclincoln;

Spring::Spring(double stiff, double length, pair<index, index> between):
  stiff_(stiff),
  length_(length),
  between_(between){

}