#include <iostream>
#include <array>
#include "separator.hpp"

int main()
{
  separator::Separator separator_solver(1.0, 1.0, 1.0);

  std::vector<Eigen::Vector3d> pointsA;
  std::vector<Eigen::Vector3d> pointsB;

  pointsA.push_back(Eigen::Vector3d(12.2, 11.4, 0));
  pointsA.push_back(Eigen::Vector3d(12, 11.6, 2));
  pointsA.push_back(Eigen::Vector3d(12.5, 11.7, 0));
  pointsA.push_back(Eigen::Vector3d(12, 11.9, 0));
  pointsA.push_back(Eigen::Vector3d(12.5, 11.4, 0));
  pointsA.push_back(Eigen::Vector3d(12.5, 11.7, 2));
  pointsA.push_back(Eigen::Vector3d(12.5, 11.4, 2));
  pointsA.push_back(Eigen::Vector3d(12.2, 11.4, 2));
  pointsA.push_back(Eigen::Vector3d(12, 11.9, 2));
  pointsA.push_back(Eigen::Vector3d(12.3, 11.8, 0));
  pointsA.push_back(Eigen::Vector3d(12, 11.6, 0));
  pointsA.push_back(Eigen::Vector3d(12.3, 11.9, 0));
  pointsA.push_back(Eigen::Vector3d(12.3, 11.9, 2));
  pointsA.push_back(Eigen::Vector3d(12.4, 11.8, 0));
  pointsA.push_back(Eigen::Vector3d(12.3, 11.8, 2));
  pointsA.push_back(Eigen::Vector3d(12.4, 11.8, 2));

  pointsB.push_back(Eigen::Vector3d(9.91, 9.86, 0.982));
  pointsB.push_back(Eigen::Vector3d(9.91, 9.86, 0.982));
  pointsB.push_back(Eigen::Vector3d(9.91, 9.86, 0.982));
  pointsB.push_back(Eigen::Vector3d(9.91, 9.86, 0.982));

  Eigen::Vector3d n;
  double d;
  bool solved = separator_solver.solveModel(n, d, pointsA, pointsB);
  std::cout << "n= " << n.transpose() << std::endl;
  std::cout << "d= " << d << std::endl;

  return 0;
};
