#include "Exporto.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>

bool Exporto::writeMatrixToOctave(Eigen::MatrixXd mat,
                                  std::filesystem::path fileName,
                                  std::string MatrixVarName, bool overwrite) {

  // Get the current time as a time_t object
  std::time_t now = std::time(0);

  // Convert time_t to local time string format
  std::string s = std::ctime(&now);
  s += ".m";

  fileName = fileName.empty() ? "ExportoOctaveFile" / std::filesystem::path(s)
                              : fileName;

  if (std::filesystem::exists(_ExportDir / fileName) && !overwrite) {
    std::cerr
        << "[ERROR] : Attempting to overwrite a file when you said not to :(\n";
    return false;
  }

  std::fstream outFile;
  outFile.open(_ExportDir / fileName, std::ios::out);

  outFile << "1;\n" << MatrixVarName << " = [";
  for (int i = 0; i < mat.rows(); i++) {
    for (int j = 0; j < mat.cols(); j++) {
      //     std::cout << "item: " << mat.row(i)(j) << "\n";

      outFile << std::setprecision(15) << mat.row(i)(j);
      if (j != mat.cols() - 1) {
        outFile << ",";
      }
    }
    //   std::cout << "next row\n";
    outFile << ";";
  }
  outFile << "];\nwhos;\n";
  std::cout << "Data written to " << _ExportDir / fileName << "\n";
  outFile.close();
  return true;
}

bool Exporto::writeMatrixToPython(Eigen::MatrixXd mat,
                                  std::filesystem::path fileName,
                                  std::string MatrixVarName, bool overwrite) {
  std::time_t now = std::time(0);
  std::string timeStr = std::ctime(&now);

  // Remove newline from timeStr
  if (!timeStr.empty() && timeStr.back() == '\n') {
    timeStr.pop_back();
  }

  fileName = fileName.empty()
                 ? std::filesystem::path("ExportoPythonFile.txt")
                 : fileName;

  // Create data filename (.txt) and script filename (.py)
  std::filesystem::path dataFileName = fileName;
  if (dataFileName.extension() != ".txt") {
    dataFileName.replace_extension(".txt");
  }

  std::filesystem::path scriptFileName = dataFileName.stem();
  scriptFileName += "_loader.py";

  // Check if files exist and overwrite flag
  if ((std::filesystem::exists(_ExportDir / dataFileName) ||
       std::filesystem::exists(_ExportDir / scriptFileName)) && !overwrite) {
    std::cerr
        << "[ERROR] : Attempting to overwrite a file when you said not to :(\n";
    return false;
  }

  // Write matrix data to text file
  std::fstream dataFile;
  dataFile.open(_ExportDir / dataFileName, std::ios::out);
  dataFile << std::setprecision(std::numeric_limits<double>::max_digits10);

  // Write matrix data row by row, space-delimited
  for (int i = 0; i < mat.rows(); i++) {
    for (int j = 0; j < mat.cols(); j++) {
      dataFile << mat(i, j);
      if (j != mat.cols() - 1) {
        dataFile << " ";
      }
    }
    dataFile << "\n";
  }
  dataFile.close();

  // Write Python loader script
  std::fstream scriptFile;
  scriptFile.open(_ExportDir / scriptFileName, std::ios::out);

  scriptFile << "#!/usr/bin/env python3\n";
  scriptFile << "# Exporto Data Export File For Python\n";
  scriptFile << "# Written on " << timeStr << "\n";
  scriptFile << "# This script loads the matrix data from " << dataFileName.filename().string() << "\n\n";
  scriptFile << "import numpy as np\n\n";
  scriptFile << "# Load matrix data from text file\n";
  scriptFile << MatrixVarName << " = np.loadtxt('" << dataFileName.filename().string() << "')\n\n";
  scriptFile << "# Display matrix information\n";
  scriptFile << "print('Matrix shape:', " << MatrixVarName << ".shape)\n";
  scriptFile << "print('Matrix dtype:', " << MatrixVarName << ".dtype)\n";
  scriptFile << "print('\\n" << MatrixVarName << ":')\n";
  scriptFile << "print(" << MatrixVarName << ")\n";

  scriptFile.close();

  std::cout << "Data written to " << _ExportDir / dataFileName << "\n";
  std::cout << "Python loader script written to " << _ExportDir / scriptFileName << "\n";

  return true;
}

bool Exporto::writePlotScript(std::filesystem::path fileName) {

  if (!std::filesystem::exists(_ExportDir / fileName)) {
    std::cerr << "[ERROR] File not found within " << _ExportDir << "\n";
    assert(std::filesystem::exists(_ExportDir / fileName));
  }

  std::fstream outFile;
  outFile.open(_ExportDir / fileName, std::ios::ate);
  outFile << "\n hold on;";
  outFile << "\n plot(arr(1,:));\n";
  outFile << "plot(arr(2,:));\n";
  outFile.close();
  return true;
}

bool Exporto::writeTimeStepVector(std::filesystem::path fileName, double T_0,
                                  double T_f, double h) {

  if (!std::filesystem::exists(_ExportDir / fileName)) {
    std::cerr << "[ERROR] File not found within " << _ExportDir << "\n";
    assert(std::filesystem::exists(_ExportDir / fileName));
  }

  std::fstream outFile;
  outFile.open(_ExportDir / fileName, std::ios::app);
  outFile << "\n"
          << fileName.stem().string() << "_"
          << "TimeSteps = " << T_0 << " : " << h << " : " << T_f << ";\n"
          << fileName.stem().string() << "_"
          << "TimeSteps =" << fileName.stem().string() << "_TimeSteps(2:end)"
          << "; \nwhos; \n ";
  outFile.close();
  return true;
}
