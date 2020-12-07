/* ----------------------------------------------------------------------------
 * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Jesus Tordesillas, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#pragma once

#include "gurobi_c++.h"
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
  GRBEnv* env_ = new GRBEnv();
  GRBModel model_ = GRBModel(*env_);

  double weight_n0_ = 0.0;
  double weight_n1_ = 0.0;
  double weight_n2_ = 0.0;

  long int num_of_LPs_run_ = 0;
  double mean_comp_time_ms_;

  double epsilon_ = 1.0;

  GRBVar n0_;  // 1st component of the normal of the plane
  GRBVar n1_;  // 2nd component of the normal of the plane
  GRBVar n2_;  // 3rd component of the normal of the plane
  GRBVar d_;   // d component of the plane
};

}  // namespace separator
