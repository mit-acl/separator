#include <iostream>
#include <array>
#include "separator.hpp"

int main()
{
  separator::Separator separator_solver;  // 0.0, 0.0, 0.0

  std::vector<Eigen::Vector3d> pointsA;
  std::vector<Eigen::Vector3d> pointsB;

  pointsA.push_back(Eigen::Vector3d(1.3901, 1.336, -5.578));
  pointsA.push_back(Eigen::Vector3d(1.7901, 1.336, -5.578));
  pointsA.push_back(Eigen::Vector3d(1.7901, 1.736, 6.422));
  pointsA.push_back(Eigen::Vector3d(1.7901, 1.336, 6.422));
  pointsA.push_back(Eigen::Vector3d(1.7901, 1.736, -5.578));
  pointsA.push_back(Eigen::Vector3d(1.3901, 1.736, 6.422));
  pointsA.push_back(Eigen::Vector3d(1.3901, 1.736, -5.578));
  pointsA.push_back(Eigen::Vector3d(1.3901, 1.336, 6.422));

  pointsB.push_back(Eigen::Vector3d(4.4376, 4.6171, 1));
  pointsB.push_back(Eigen::Vector3d(4.4376, 4.6171, 1));
  pointsB.push_back(Eigen::Vector3d(4.4376, 4.6171, 1));
  pointsB.push_back(Eigen::Vector3d(4.4376, 4.6171, 1));

  Eigen::Vector3d n;
  double d;
  bool solved = separator_solver.solveModel(n, d, pointsA, pointsB);
  std::cout << "solved= " << solved << std::endl;
  std::cout << "n= " << n.transpose() << std::endl;
  std::cout << "d= " << d << std::endl;

  return 0;
};
