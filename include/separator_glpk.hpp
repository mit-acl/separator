/* ----------------------------------------------------------------------------
 * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Jesus Tordesillas, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

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
  Separator();

  bool solveModel(Eigen::Vector3d& solution, double& d, const std::vector<Eigen::Vector3d>& pointsA,
                  const std::vector<Eigen::Vector3d>& pointsB);

  bool solveModel(Eigen::Vector3d& solutionN, double& solutionD,
                  const Eigen::Matrix<double, 3, Eigen::Dynamic>& pointsA,
                  const Eigen::Matrix<double, 3, Eigen::Dynamic>& pointsB);

  long int getNumOfLPsRun();

  double meanSolveTimeMs();

private:
  glp_prob* lp_;
  int ia_[10000], ja_[10000];  // TODO
  double ar_[10000];           // TODO
  double weight_n0_ = 0.0;
  double weight_n1_ = 0.0;
  double weight_n2_ = 0.0;

  glp_smcp params_;

  long int num_of_LPs_run_ = 0;
  double mean_comp_time_ms_;
};

}  // namespace separator
