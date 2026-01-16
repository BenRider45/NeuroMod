#include <Eigen/Dense>
#include <filesystem>
#include <string>
class Exporto {
public:
  Exporto(std::filesystem::path ExportDir) : _ExportDir(ExportDir) {
    assert(std::filesystem::exists(ExportDir));
  };

  bool writeMatrixToOctave(Eigen::MatrixXd mat, std::filesystem::path fileName,
                           std::string MatrixVarName, bool overwrite = false);

  bool writeMatrixToPython(Eigen::MatrixXd mat, std::filesystem::path fileName,
                           std::string MatrixVarName, bool overwrite = false);

  bool writeTimeStepVector(std::filesystem::path fileName, double T_0,
                           double T_f, double h);
  bool writePlotScript(std::filesystem::path fileName);

private:
  std::filesystem::path _ExportDir;
};
