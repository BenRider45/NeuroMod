
#include "HVCRA.hpp"
#include <iostream>
#include <limits>

Eigen::MatrixXd HVCRA::SimDataToMatrix(SimData data) {
  Eigen::MatrixXd m(data.front().rows(), data.size());
  std::cout << "m: " << m.rows() << "X" << m.cols() << "\n";
  int i = 0;
  for (const auto col : data) {
    m.col(i) = col;
    i++;
  }

  return m;
}

// double I_sL(double V_s) { return _consts.G_L * (V_s - _consts.E_L); }
namespace {
double m_alpha(double V_s) {
  return (.1 * (V_s + 40)) / (1 - std::exp(-.1 * (V_s + 40)));
}

double m_beta(double V_s) { return 4 * std::exp(-0.0556 * (V_s + 65)); }

double tao_m(double V_s) { return 1 / (m_alpha(V_s) + m_beta(V_s)); }

double inf_m(double V_s) {
  return (m_alpha(V_s)) / (m_alpha(V_s) + m_beta(V_s));
}

} // namespace

Eigen::MatrixXd HVCRA::GimmeCurrents(Eigen::MatrixXd Currents,
                                     Eigen::MatrixXd GatingVars,
                                     Eigen::MatrixXd CalciumConc) {
  // Returns Is_l, I_sNa, I_sK, I_dL, I_dCa, I_dCaK
  //
  // Gating Vars of form h,n,r,c
  //
  std::cout << Currents.cols() << "==" << GatingVars.cols()
            << "==" << CalciumConc.cols() << "\n";
  assert(Currents.cols() == GatingVars.cols() &&
         GatingVars.cols() == CalciumConc.cols());

  Eigen::MatrixXd output(6, Currents.cols());

  for (int i = 0; i < Currents.cols(); i++) {
    double V_s = Currents(0, i);
    double V_d = Currents(1, i);
    double h = GatingVars(0, i);
    double n = GatingVars(1, i);
    double r = GatingVars(2, i);
    double c = GatingVars(3, i);
    double Ca = CalciumConc(0, i);

    output(0, i) = -_consts.G_L * (V_s - _consts.E_L);

    output(1, i) = -_consts.G_s_Na * m_inf(V_s) * m_inf(V_s) * m_inf(V_s) *
                   GatingVars(0, i) * (V_s - _consts.E_Na);

    output(2, i) = -_consts.G_Kdr * n * n * n * n * (V_s - _consts.E_K);
    output(3, i) = -_consts.G_L * (V_d - _consts.E_L);
    output(4, i) = -_consts.G_Ca * r * r * (V_d - _consts.E_Ca);
    output(5, i) =
        -((_consts.G_CaK * c) / (1.0 + (6.0 / Ca))) * (V_d - _consts.E_Ca);
  }
  return output;
}

State HVCRA::f(double t, State &y) {

  // y(0) V_
  // y(4) r
  // y(5) c
  // y(6) Ca
  // y(7) g_s_exc
  // y(8) g_d_exc
  // y(9) g_s_inh
  // y(10) g_d_inh
  // y(11) m

  double I_sL, I_sNa, I_sKdr, I_sExc, I_sInh;
  I_sL = -_consts.G_L * (y(0) - _consts.E_L);
  double m_inf_val = m_inf(y(0));
  I_sNa = -_consts.G_s_Na * m_inf_val * m_inf_val * m_inf_val * y(2) *
          (y(0) - _consts.E_Na);

  I_sKdr = -_consts.G_Kdr * y(3) * y(3) * y(3) * y(3) * (y(0) - _consts.E_K);

  I_sExc = -y(7) * y(0);
  I_sInh = -y(9) * (y(0) - _consts.E_I);

  double I_Ts = I_sL + I_sNa + I_sKdr;
  double I_dL, I_dCa, I_dCaK, I_dExc, I_dInh;

  I_dL = -_consts.G_L * (y(1) - _consts.E_L);
  I_dCa = -_consts.G_Ca * y(4) * y(4) * (y(1) - _consts.E_Ca);

  I_dCaK =
      (-(_consts.G_CaK * y(5)) / (1.0 + (6.0 / y(6)))) * (y(1) - _consts.E_K);

  I_dExc = -y(8) * y(1);
  I_dInh = -y(10) * (y(1) - _consts.E_I);
  double I_Td = I_dL + I_dCa + I_dCaK + I_dExc + I_dInh;
  if (t < 40) {

    //   std::cout << "_I_sExt(" << t << "): " << _I_sExt(t) << "\n";
    //   std::cout << "_I_dExt: " << _I_dExt(t) << "\n";
  }
  State out(11);
  out(0) =
      (I_Ts / _consts.C_m +
       10e5 *
           (

               (_I_sExt(t) / (_consts.C_m * _consts.A_s))

               + ((y(1) - y(0)) / (_consts.R_c * _consts.C_m * _consts.A_s))));

  out(1) =
      I_Td / _consts.C_m +

      10e5 * (

                 (_I_dExt(t) / (_consts.C_m * _consts.A_d)) +

                 ((y(0) - y(1)) / (_consts.R_c * _consts.C_m * _consts.A_d))

             );

  out(2) = (h_inf(y(0)) - y(2)) / tao_h(y(0));

  out(3) = (n_inf(y(0)) - y(3)) / tao_n(y(0));

  out(4) = (r_inf(y(1)) - y(4)) / tao_r(y(1));

  out(5) = (c_inf(y(1)) - y(5)) / tao_c(y(1));

  out(6) = .1 * I_dCa - .02 * y(6);

  out(7) = -(1.0 / _consts.tao) * y(7);

  out(8) = -(1.0 / _consts.tao) * y(8);

  out(9) = -(1.0 / _consts.tao) * y(9);

  out(10) = -(1.0 / _consts.tao) * y(10);

  // out(11) = m_alpha() * (1 - y(11)) - m_beta(
  for (const auto num : out) {
    if (num > std::numeric_limits<double>::max()) {
      std::cerr << "ERROR: Double Max Reached\n";
      assert(false);
    }
  }
  return out;
}

bool HVCRA::verify_state(const State &y) { return y.size() == ode_num; }

void HVCRA::post_step_process(const double t, State &y) {
  if (y(0) > -10e-3) {
    // std::cout << "Somatic Spike at t=" << t << "\n";
  }

  if (y(1) > -10e-3) {
    // std::cout << "Dendritic Spike at t=" << t << "\n";
  }
}
