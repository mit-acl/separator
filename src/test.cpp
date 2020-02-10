#include <iostream>
#include <array>
#include "separator.hpp"

int main()
{
  separator::Separator separator_solver(0.0, 0.0, 0.0);

  std::vector<Eigen::Vector3d> pointsA;
  std::vector<Eigen::Vector3d> pointsB;

  pointsA.push_back(Eigen::Vector3d(1, 0, 0));
  pointsA.push_back(Eigen::Vector3d(0.5, 0, 0));
  pointsB.push_back(Eigen::Vector3d(-1, 0, 0));

  Eigen::Vector3d n;
  bool solved = separator_solver.solveModel(n, pointsA, pointsB);
  std::cout << "n= " << n.transpose() << std::endl;

  return 0;
};
