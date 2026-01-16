#include "Rk4.hpp"

Eigen::MatrixXd RK4::Simulate(OdeAble &system, State y_0, double t_0,
                              double t_f) {
  // int timesteps = static_cast<int>(std::floor((t_f - t_0) / _h)) + 1;
  int timesteps = (t_f - t_0) / _h;
  assert(system.verify_state(y_0));

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> data(system.ode_num,
                                                             timesteps);

  data.col(0) = y_0;

  // int n = static_cast<int>(std::floor((t_f - t_0) / _h)) + 1;
  Eigen::VectorXd timeSteps = Eigen::VectorXd::LinSpaced(timesteps, t_0, t_f);
  double t = timeSteps(1);

  for (int i = 0; i < timeSteps.size() - 1; i++) {

    t = timeSteps(i);
    State last = data.col(i);
    State y_i = Step(system, t, last);
    system.post_step_process(t, y_i);
    data.col(i + 1) = y_i;
  }
  data.conservativeResize(system.ode_num + 1, Eigen::NoChange_t());
  data.row(system.ode_num) = timeSteps;
  return data;
}

State RK4::Step(OdeAble &system, double t, State &y) {
  // RK4 method here!
  //% Calculate the four stages
  // K1 = feval(f, t(i),y(:,i) );
  // K2 = feval(f, t(i)+.5*h, y(:,i)+.5*h*K1);
  // K3 = feval(f, t(i)+.5*h, y(:,i)+.5*h*K2);
  // K4 = feval(f, t(i)+h, y(:,i)+h*K3 );
  //% Evaluate approximate solution at next step
  // y(:,i+1) = y(:,i) + h/6 *(K1+2*K2+2*K3+K4);
  // end
  State K1 = system.f(t, y);
  State K2_y = y + .5 * _h * K1;

  State K2 = system.f(t + .5 * _h, K2_y);
  State K3_y = y + .5 * _h * K2;
  State K3 = system.f(t + .5 * _h, K3_y);
  //
  State K4_y = y + _h * K3;
  State K4 = system.f(t + _h, K4_y);
  //

  //  State Y_1 = y;
  //  State Y_2 = y + _h / 2 * system.f(t, Y_1);
  //  State Y_3 = y + _h / 2 * system.f(t + (_h / 2), Y_2);
  //  State Y_4 = y + _h * system.f(t + (_h / 2), Y_3);
  //
  //  return y +
  //         (_h / 6) * (system.f(t, Y_1) + 2 * system.f(t + (_h / 2), Y_2) +
  //                     2 * system.f(t + (_h / 2), Y_3) + system.f(t + _h,
  //                     Y_4));
  return y + ((_h / 6) * (K1 + 2 * K2 + 2 * K3 + K4));
}
