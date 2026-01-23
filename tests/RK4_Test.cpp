#include "BCNNeuroModel/OdeAble.hpp"
#include "BCNNeuroModel/Rk4.hpp"
#include <cmath>
#include <gtest/gtest.h>

class EasyOde : public OdeAble {
public:
  EasyOde(double ode_num) : OdeAble(1) {};

  State f(double t, State &y) override {
    State output(1);

    output << -y * y;

    return output;
  }

  bool verify_state(State const &y) override { return true; }

  void post_step_process(double t, State &y) override { return; }

  // yeah wirk krirkenueninely yessirthanks mr code':)
};
double getRate(double e_h, double e_2h) { return std::log2(e_2h / e_h); }
double getActual(double t) { return 1 / t; }
TEST(RK4, isRate4) {
  ASSERT_TRUE(true);

  // Implementing Verification method from Ch16 Ascher, Greif. (A First Course
  // In Numerical Methods)
  double h_init = 1e-3;
  std::vector<double> steps = {.002, .005, .01, .02, .05, .1, .2};
  double t_f = 10.0;
  double t_0 = 1.0;
  EasyOde easy(1.0);
  State easy_init(1);
  easy_init << 1;

  for (double h : steps) {
    RK4 solver(h);
    RK4 solverDouble(h * 2);

    Eigen::MatrixXd output = solver.Simulate(easy, easy_init, t_0, t_f);
    Eigen::MatrixXd outputDoubled =
        solverDouble.Simulate(easy, easy_init, 1.0, 10.0);

    std::cerr << "Output1: " << output.rows() << "," << output.cols() << "\n";

    std::cerr << "Output1: " << output.col(output.cols() - 1) << "\n";
    std::cerr << "Output2: " << outputDoubled.col(outputDoubled.cols() - 1)
              << "\n"
              << std::flush;
    double Errh = std::abs(output(0, output.cols() - 1) - getActual(t_f));
    std::cerr << "Error for timestep " << h << ": " << Errh << "\n";

    double Err2h =
        std::abs(outputDoubled(0, outputDoubled.cols() - 1) - getActual(t_f));
    std::cerr << "Error for timestep " << h * 2 << ": " << Err2h << "\n";
    std::cerr << "Rate: " << getRate(Errh, Err2h);
    EXPECT_TRUE(getRate(Errh, Err2h) >= 3.8);
  }
};
