#include <Eigen/Dense>
#include <functional>

namespace CustomTypes {

typedef Eigen::VectorXd EigVector;
typedef std::function<EigVector(double, const VectorXd &)> OdeFunc;
typedef std::function<double(double)> ExtCurrent;

} // namespace CustomTypes
