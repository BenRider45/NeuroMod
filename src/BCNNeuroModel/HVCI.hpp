#pragma once

#include "OdeAble.hpp"
#include <Eigen/Dense>
#include <cmath>
#include <functional>
struct HVCI_CONSTANTS {
  HVCI_CONSTANTS(double E_l, double E_Na, double E_K, double E_I, double E_E,
                 double G_L, double G_Na, double G_Kdr, double G_KHT, double A,
                 double C_m, double tau_exc, double tau_inh)
      : _E_l(E_l), _E_Na(E_Na), _E_K(E_K), _E_I(E_I), _E_E(E_E), _G_L(G_L),
        _G_Na(G_Na), _G_Kdr(G_Kdr), _G_KHT(G_KHT), _A(A), _C_m(C_m),
        _tau_exc(tau_exc), _tau_inh(tau_inh) {};
  double _E_l;
  double _E_Na;
  double _E_K;
  double _E_I;
  double _E_E;
  double _G_L;
  double _G_Na;
  double _G_Kdr;
  double _G_KHT;
  double _A;
  double _C_m;
  double _tau_exc;
  double _tau_inh;
};

constexpr int HVCI_NUM_ODE = 7;

class HVCI : public OdeAble {
public:
  HVCI(HVCI_CONSTANTS consts, std::function<double(double)> I_ext, double G_E_0,
       double G_I_0, double id)
      : OdeAble(HVCI_NUM_ODE), _consts(consts), _I_ext(I_ext), id(id) {}

  State f(double t, State &y) override;

  bool verify_state(State const &y) override;

  void post_step_process(const double t, State &y) override;

  int id;

private:
  std::function<double(double)> _I_ext;

  double a_n(double V) {
    return (0.15 * (V + 15)) / (1 - std::exp(-(V + 15) / 10));
  };

  HVCI_CONSTANTS _consts;

  double a_m(double V) { return (V + 22) / (1 - std::exp(-(V + 22) / 10)); };

  double a_h(double V) { return (1) / (1 + std::exp(-V / 5)); }

  double w_inf(double V) { return 1 / (1 + std::exp(-V / 5)); };

  double tau_w(double V) { return 1.0; };

  double b_n(double V) { return .2 * std::exp(-(V + 25) / 80); };

  double b_m(double V) { return 40 * std::exp(-(V + 47) / 18); };

  double b_h(double V) { return 10 / (1 + std::exp(-(V + 4) / 10)); };

  double I_L(double V);

  double I_Kdr(double V, double n);

  double I_KHT(double V, double w);

  double I_Na(double V, double m, double h);

  double I_exec(double V, double G_E);

  double I_inh(double V, double G_I);
};
