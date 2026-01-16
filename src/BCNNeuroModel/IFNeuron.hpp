#pragma once
#include "OdeAble.hpp"
#include <cmath>
#include <functional>
struct IFNEURON_CONSTANTS {
  IFNEURON_CONSTANTS(double I_0, double L, double Tau)
      : _I_0(I_0), _L(L), _Tau(Tau) {};

  double _I_0;
  double _L;
  double _Tau;
};
class IFNeuron : public OdeAble {
public:
  IFNeuron(IFNEURON_CONSTANTS consts, std::function<double(double)> Ifunc)
      : OdeAble(1), _consts(consts), I(Ifunc) {}
  virtual State f(double t, State &y) override {
    // y(0) = V

    State out(1);

    out(0) = (_consts._L - y(0) + _consts._I_0) / _consts._Tau;
    return out;
  };

  virtual void post_step_process(const double t, State &y_i) override {

    return;
  };

  virtual bool verify_state(const State &y_i) override {
    return y_i.size() == ode_num;
  };

private:
  IFNEURON_CONSTANTS _consts;
  std::function<double(double)> I;
};
