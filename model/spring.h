#ifndef _SPRING_H_
#define _SPRING_H_

#include <utility>

namespace zyclincoln{

  struct Spring{
  public:
    Spring(const double stiff, const double length, const std::pair<unsigned int, unsigned int> between);
    double stiff_;
    double length_;
    std::pair<unsigned int, unsigned int> between_;
  };

}

#endif