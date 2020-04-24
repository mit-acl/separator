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

  // void deleteModel();

private:
  glp_prob* lp_;
  int ia_[1 + 1000], ja_[1 + 1000];
  double ar_[1 + 1000];
  double weight_n1_ = 1.0;
  double weight_n2_ = 1.0;
  double weight_n3_ = 1.0;

  glp_smcp params_;
};

}  // namespace separator
