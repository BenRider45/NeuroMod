#include "HVCI.hpp"

double HVCI::I_L(double V) { return -_consts._G_L * (V - _consts._E_l); }

double HVCI::I_Kdr(double V, double n) {
  return -_consts._G_Kdr * n * n * n * n * (V - _consts._E_K);
}

double HVCI::I_KHT(double V, double w) {
  return -_consts._G_KHT * w * (V - _consts._E_K);
}

double HVCI::I_Na(double V, double m, double h) {
  return -_consts._G_Na * m * m * m * h * (V - _consts._E_Na);
}

double HVCI::I_exec(double V, double G_E) { return -G_E * (V - _consts._E_E); }

double HVCI::I_inh(double V, double G_I) { return -G_I * (V - _consts._E_I); }

State HVCI::f(double t, State &y) {
  // y(0) V
  // y(1) n
  // y(2) m
  // y(3) h
  // y(4) w
  // y(5) g_exc
  // y(6) g_inh

  State output(this->ode_num);

  double I_t, I_L, I_Kdr, I_KHT, I_Na, I_exc, I_inh;

  I_L = HVCI::I_L(y(0));

  I_Kdr = HVCI::I_Kdr(y(0), y(1));

  I_KHT = HVCI::I_KHT(y(0), y(4));

  I_Na = HVCI::I_Na(y(0), y(2), y(3));

  I_exc = HVCI::I_exec(y(0), y(5));

  I_inh = HVCI::I_inh(y(0), y(6));

  I_t = I_L + I_Kdr + I_KHT + I_Na + I_exc + I_inh;

  output(0) =
      (I_t / _consts._C_m) + 10e4 * (_I_ext(t) / (_consts._C_m * _consts._A));

  output(1) = ((a_n(y(0)) * (1 - y(1))) - b_n(y(0)) * y(1));

  output(2) = ((a_m(y(0)) * (1 - y(2))) - b_m(y(0)) * y(2));

  output(3) = ((a_h(y(0)) * (1 - y(3))) - b_h(y(0)) * y(3));

  output(4) = (w_inf(y(0)) - y(4)) / tau_w(y(0));

  output(5) = -y(5) / _consts._tau_exc;

  output(6) = -y(6) / _consts._tau_inh;

  return output;
}

bool HVCI::verify_state(State const &y) { return y.size() == this->ode_num; }

void HVCI::post_step_process(const double t, State &y) { return; }
