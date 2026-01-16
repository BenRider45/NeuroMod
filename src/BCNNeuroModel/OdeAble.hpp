#pragma once

#include <Eigen/Dense>
#include <cassert>
typedef Eigen::VectorXd State;

class OdeAble {
public:
  virtual ~OdeAble() {};

  OdeAble(int ode_num) : ode_num(ode_num) { assert(ode_num > 0); };
  virtual State f(double t, State &y) = 0;
  virtual bool verify_state(State const &y) = 0;
  virtual void post_step_process(double t, State &y) = 0;
  int ode_num;
};
