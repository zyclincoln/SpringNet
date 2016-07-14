#ifndef _IMPLICIT_EULER_INTEGRATOR_H_
#define _IMPLICIT_EULER_INTEGRATOR_H_

#include "integrator.h"

namespace zyclincoln{

  class ImplicitEulerIntegrator : public Integrator{
  public:
    virtual void next_frame(AbstractSystem& system);
  };

}

#endif