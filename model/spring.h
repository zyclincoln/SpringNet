#ifndef _SPRING_H_
#define _SPRING_H_

namespace zyclincoln{

  struct Spring{
  public:
    Spring(double stiff, double length, std::pair<unsigned int, unsigned int> between);
    double stiff_;
    double length_;
    std::pair<unsigned int, unsigned int> between_;
  };

}

#endif