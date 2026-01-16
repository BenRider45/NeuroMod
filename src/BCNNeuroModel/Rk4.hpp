#pragma once
#include "OdeAble.hpp"
#include <Eigen/Dense>
#include <cassert>
#include <iostream>
class RK4 {
public:
  RK4(double h) : _h(h) {
    if (h <= 0.0) {
      std::cerr << "[ERROR] time step must be greater than 0 :(\n";
      assert(h > 0);
    }
  };
  Eigen::MatrixXd Simulate(OdeAble &system, State y_0, double t_0, double t_f);

private:
  double _h;

  State Step(OdeAble &system, double t, State &y);
};
