#pragma once

#include <glpk.h> /* GNU GLPK linear/mixed integer solver */
#include <array>
#include <vector>
#include <Eigen/Dense>

namespace separator
{
class Separator
{
public:
  Separator();  // double weight_n1, double weight_n2, double weight_n3

  bool solveModel(Eigen::Vector3d& solution, double& d, const std::vector<Eigen::Vector3d>& pointsA,
                  const std::vector<Eigen::Vector3d>& pointsB);

  bool solveModel(Eigen::Vector3d& solutionN, double& solutionD,
                  const Eigen::Matrix<double, 3, Eigen::Dynamic>& pointsA,
                  const Eigen::Matrix<double, 3, Eigen::Dynamic>& pointsB);

  // void deleteModel();

private:
  glp_prob* lp_;
  int ia_[10000], ja_[10000];  // TODO
  double ar_[10000];           // TODO
  double weight_n1_ = 1.0;
  double weight_n2_ = 1.0;
  double weight_n3_ = 1.0;

  glp_smcp params_;
};

}  // namespace separator
