// #include "CustomTypes.hpp"
#pragma once
#include "OdeAble.hpp"
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <functional>

// TODO This namespace thing kinda sucks, undo that
static double TIME_STEP = 1e-3;

struct HVCRA_CONSTANTS {
  HVCRA_CONSTANTS(double C_m, double A_s, double A_d, double R_c, double G_L,
                  double G_dL, double G_Ca, double G_s_Na, double E_Ca,
                  double E_L, double E_Na, double G_Kdr, double E_K, double E_I,
                  double tao, double G_CaK, double spike_conductivity_G) {
    this->C_m = C_m;
    this->A_s = A_s;
    this->A_d = A_d;
    this->R_c = R_c;
    this->G_L = G_L;
    this->G_dL = G_dL;
    this->G_Ca = G_Ca;
    this->G_s_Na = G_s_Na;
    this->E_Ca = E_Ca;
    this->E_L = E_L;
    this->E_Na = E_Na;
    this->G_Kdr = G_Kdr;
    this->E_K = E_K;
    this->E_I = E_I;
    this->tao = tao;
    this->G_CaK = G_CaK;
    this->spike_conductivity_G = spike_conductivity_G;
  };

  double C_m;    // Membrane Capacitance
  double A_s;    // Somatic Compartment Area
  double A_d;    // Dendritic Compartment Area
  double R_c;    // Coupling Resistance btwn compartments
  double G_L;    // Leak Conductance
  double G_dL;   // Leak Conductance (Dendritic Component)
  double G_Ca;   // Ca^++ current Conductance
  double G_s_Na; // Na^+ currrent conductance
  double E_Ca;   // CA^++ high threshold current reversal potential
  double E_L;    // Leak current reversal potential
  double E_Na;   // Na^+ current reversal potential
  double G_Kdr;
  double G_CaK; // Delay-Rectified K+ current conductance
  double E_K;   // Delay-Rectified K+ current reversal potential
  double E_I;   // Inhibitory Synaptic current Reversal Potential
  double tao;   // Synaptic Time Constant
  double spike_conductivity_G;
};

static int NUM_ODE = 11;

typedef std::function<State(double, const State &)> OdeFunc;
typedef std::vector<State> SimData;

class HVCRA : public OdeAble {

public:
  HVCRA(HVCRA_CONSTANTS consts, std::function<double(double)> I_sExt,
        std::function<double(double)> I_dExt, double G, int id)
      : OdeAble(NUM_ODE), _consts(consts), id(id), _I_sExt(I_sExt),
        _I_dExt(I_dExt) {}

  State f(double t, State &y) override;

  bool verify_state(State const &y) override;

  void post_step_process(const double t, State &y) override;

  // y(0) V_s
  // y(1) V_d
  // y(2) h
  // y(3) n
  // y(4) r
  // y(5) c
  // y(6) Ca
  // y(7) g_s_exc
  // y(8) g_d_exc
  // y(9) g_s_inh
  // y(10) g_d_inh
  SimData simulate(double t_0, double t_f, double h, State y_0);
  Eigen::MatrixXd SimDataToMatrix(SimData data);
  int id;

private:
  State step(double &t, double h, State y);
  HVCRA_CONSTANTS _consts;
  std::function<double(double)> _I_sExt;
  std::function<double(double)> _I_dExt;

  double tao_h(double V_s) {
    return 0.1 + 0.75 / (1.0 + std::exp((V_s + 40.5) / 6.0));
  };

  double h_inf(double V_s) { return 1.0 / (1 + std::exp((V_s + 45.0) / 7.0)); };

  double tao_n(double V_s) {
    return 0.1 + 0.5 / (1.0 + std::exp((V_s + 27.0) / 15.0));
  };
  double n_inf(double V_s) {
    return 1.0 / (1.0 + std::exp(-(V_s + 35.0) / 10.0));
  };

  double tao_r(double V_d) { return 1.0; };
  double r_inf(double V_d) {
    return 1.0 / (1.0 + std::exp(-(V_d + 5.0) / 10.0));
  };

  double tao_c(double V_d) { return 10.0; };
  double c_inf(double V_d) {
    return 1.0 / (1.0 + std::exp(-(V_d - 10.0) / 7.0));
  };

  double m_inf(double V_s) {
    return 1.0 / (1.0 + std::exp(-(V_s + 30.0) / 9.5));
  };
};

// namespace HVCRA
