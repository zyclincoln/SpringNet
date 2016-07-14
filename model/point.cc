#include "point.h"

using namespace zyclincoln;
using namespace Eigen;

Point::Point(const Vector3d &position, const Vector3d &speed, const double mass):
  position_(position),
  speed_(speed),
  mass_(mass){

}