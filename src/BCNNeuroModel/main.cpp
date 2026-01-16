#include <iomanip>
#include <iostream>

#include "Exporto.hpp"
#include "HVCRA.hpp"
#include "IFNeuron.hpp"
#include "Rk4.hpp"
double SomaExternalCurrent(double t) { return t <= 20e-3 ?: 0; }

double DendriteExternalCurrent(double t) {
  return t >= 20.0 && t <= 40.0 ? 1.6e9 : 0;
}

int gen_HVCRA() {

  HVCRA_CONSTANTS consts(1, 5000, 10000, 5.5, .1, 0.1, 55, 60, 120, -80, 55, 8,
                         -90, -80, 5, 150, -10);

  HVCRA neuron(consts, SomaExternalCurrent, DendriteExternalCurrent, .01, -1);

  double t_0 = 0, t_f = 100, h = .001;

  std::cout << "Neuron id: " << neuron.id << "\n";
  State y_0(11);

  y_0 << -79.9887, -79.9660, 0.9933, 0.0110, 0.0006, 0.0000, 0.0152, 0, 0, 0, 0;
  // y_0 << -75, -75, 0.0, 0.0, 0.0, 0.0000, 0.0, 0, 0, 0, 0;
  //    -75,     // V_s: Somatic voltage (resting potential ~-65 mV)
  //        -75, // V_d: Dendritic voltage (resting potential ~-65 mV)
  //        0.0, // h: Sodium inactivation (near fully inactivated at rest)
  //        0.0, // n: Potassium activation (minimally activated at rest)
  //        0.0, // r: Calcium activation (closed at rest)
  //        0.0, // c: Ca-K activation (closed at rest)
  //        0.0, // Ca: Resting calcium concentration (~50 nM)
  //        0.0, // g_s_exc: No initial excitatory input to soma
  //        0.0, // g_d_exc: No initial excitatory input to dendrite
  //        0.0, // g_s_inh: No initial inhibitory input to soma
  //        0.0; // g_d_inh: No initial inhibitory input to dendrite
  //
  //
  //
  RK4 solver(h);
  Eigen::MatrixXd data = solver.Simulate(neuron, y_0, t_0, t_f);
  //  for (auto vec : data) {
  //    std::cout << "Vec: " << vec << "\n";
  //  }
  // Eigen::MatrixXd m = neuron.SimDataToMatrix(data);
  // std::cout << "MatrixXd: " << m << "\n";
  std::filesystem::path p = "./";
  Exporto exp(p);
  // Eigen::MatrixXd m1(3, 3);
  // m1 << 1, 2, 3, 4, 5, 6, 7, 8, 9;
  exp.writeMatrixToPython(data, "HVCRA_DATA.py", "HVCRA_DATA", true);
  // exp.writePlotScript("Balls.m");

  return 0;
}

double x_squared(double x) { return x * x; }

double IF_Ifunc(double t) { return 0; }

int IFNeuronStuff() {
  std::vector<double> timesteps = {.01, .02, .04, .1, .2, .001};
  IFNEURON_CONSTANTS ifConsts(5, -70, 20);
  IFNeuron neuron(ifConsts, IF_Ifunc);
  State y_0_if(1);
  y_0_if(0) = -70;
  int i = 0;
  std::vector<double> endVals;

  for (double h : timesteps) {
    RK4 solver(h);
    Eigen::MatrixXd data = solver.Simulate(neuron, y_0_if, 0, 40);
    endVals.push_back(data(0, data.cols() - 1));
    std::cout << "data: " << data.rows() << "," << data.cols() << "\n";
    Exporto exp("./IFNData");
    std::string varName = "IFNDataH" + std::to_string(i);
    std::string fileName = varName + ".m";
    exp.writeMatrixToPython(data, fileName, varName, true);
    i++;
  }

  double goodVal = endVals.back();

  std::vector<double> errVals;

  std::cout << "[";
  for (int i = 0; i < endVals.size() - 1; i++) {
    errVals.push_back(std::abs(endVals[i] - goodVal));
    std::cout << timesteps[i] << ",";
  }
  std::cout << "]\n";

  std::cout << "[";

  for (double v : errVals) {
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10)
              << v << ",";
  }

  std::cout << "]";
  //  RK4 solver(1e-5);
  //  Eigen::MatrixXd data = solver.Simulate(neuron, y_0_if, 0, 100e-3);
  //  Exporto exp("./IFNData");
  //  exp.writeMatrixToOctave(data, "IFNDataH0.m", "IFNDataH0", true);
  //  exp.writeTimeStepVector("IFNDataH0.m", 0.0, 100e-3, 1e-5);
  //
  //  RK4 solver1(1e-4);
  //  data = solver1.Simulate(neuron, y_0_if, 0, 100e-3);
  //  exp.writeMatrixToOctave(data, "IFNDataH1.m", "IFNDataH1", true);
  //  exp.writeTimeStepVector("IFNDataH1.m", 0.0, 100e-3, 1e-4);
  //
  //  RK4 solver2(1e-3);
  //  data = solver2.Simulate(neuron, y_0_if, 0, 100e-3);
  //  exp.writeMatrixToOctave(data, "IFNDataH2.m", "IFNDataH2", true);
  //  exp.writeTimeStepVector("IFNDataH2.m", 0.0, 100e-3, 1e-3);
  //
  //  RK4 solver3(1e-2);
  //  data = solver3.Simulate(neuron, y_0_if, 0, 100e-3);
  //  exp.writeMatrixToOctave(data, "IFNDataH3.m", "IFNDataH3", true);
  //  exp.writeTimeStepVector("IFNDataH3.m", 0.0, 100e-3, 1e-2);
  return 0;
}

int main(int argc, char *argv[]) {
  //  double h = 1;
  //  double t_0 = 0.0;
  //  double t_1 = 5;
  //
  //  double t = t_0;
  //  int timeSteps = std::ceil((t_1 - t_0) / h);
  //  Eigen::MatrixXd vals(1, timeSteps);
  //
  //  for (int i = 0; i < ((t_1 - t_0) / h); i++) {
  //    vals(i) = x_squared(t);
  //    t += h;
  //  }
  //
  //  Exporto exp("./FuncTest");
  //
  //  exp.writeMatrixToPython(vals, "XSqu.py", "data", true);
  //
  // IFNeuronStuff();
  gen_HVCRA();
  return 0;
}
