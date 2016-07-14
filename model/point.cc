#include <Eigen/Dense>
#include "point.h"

using namespace zyclincoln;

Point::Point(Vector3d &position, Vector3d &speed, double mass):
  position_(position),
  speed_(speed),
  mass_(mass){

}