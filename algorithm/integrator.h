#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include <Eigen/Dense>

#include "../model/abstract_system.h"

namespace zyclincoln{

  class Integrator{
  public:
    virtual void next_frame(AbstractSystem& system) = 0;
  }

}

#endif