/* ----------------------------------------------------------------------------
 * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Jesus Tordesillas, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include <time.h>
#include <math.h>

#include <chrono>

#include "separator_gurobi.hpp"
#include <iostream>

namespace separator
{
Separator::Separator()  // double weight_n1, double weight_n2, double weight_n3
{
  n0_ = model_.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS);  //, coeff[i] + std::to_string(t)
  n1_ = model_.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS);  //, coeff[i] + std::to_string(t)
  n2_ = model_.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS);  //, coeff[i] + std::to_string(t)
  d_ = model_.addVar(-GRB_INFINITY, GRB_INFINITY, 0, GRB_CONTINUOUS);   //, coeff[i] + std::to_string(t)

  model_.set("OutputFlag", std::to_string(0));  // 1 if you want verbose, 0 if not

  epsilon_ = 1.0;
};

long int Separator::getNumOfLPsRun()
{
  return num_of_LPs_run_;
}

double Separator::meanSolveTimeMs()
{
  return mean_comp_time_ms_;
}

bool Separator::solveModel(Eigen::Vector3d& solutionN, double& solutionD, const std::vector<Eigen::Vector3d>& pointsA,
                           const std::vector<Eigen::Vector3d>& pointsB)
{
  Eigen::Matrix<double, 3, Eigen::Dynamic> pointsA_matrix(3, pointsA.size());
  for (int i = 0; i < pointsA.size(); i++)
  {
    pointsA_matrix.col(i) = pointsA[i];
  }

  Eigen::Matrix<double, 3, Eigen::Dynamic> pointsB_matrix(3, pointsB.size());

  for (int i = 0; i < pointsB.size(); i++)
  {
    pointsB_matrix.col(i) = pointsB[i];
  }

  return solveModel(solutionN, solutionD, pointsA_matrix, pointsB_matrix);
}

bool Separator::solveModel(Eigen::Vector3d& solutionN, double& solutionD,
                           const Eigen::Matrix<double, 3, Eigen::Dynamic>& pointsA,
                           const Eigen::Matrix<double, 3, Eigen::Dynamic>& pointsB)
{
  auto start_time = std::chrono::high_resolution_clock::now();

  // std::cout << "pointsA_matrix.cols()=" << pointsA.cols() << std::endl;
  // std::cout << "pointsB_matrix.cols()=" << pointsB.cols() << std::endl;

  ////////////////////////////  RESET (except for the variables)
  ////////////////////////////

  GRBConstr* c = 0;
  c = model_.getConstrs();
  for (int i = 0; i < model_.get(GRB_IntAttr_NumConstrs); ++i)
  {
    model_.remove(c[i]);
  }

  GRBQConstr* cq = 0;
  cq = model_.getQConstrs();
  for (int i = 0; i < model_.get(GRB_IntAttr_NumQConstrs); ++i)
  {
    model_.remove(cq[i]);
  }

  GRBGenConstr* gc = 0;
  gc = model_.getGenConstrs();
  for (int i = 0; i < model_.get(GRB_IntAttr_NumGenConstrs); ++i)
  {
    model_.remove(gc[i]);
  }

  // GRBVar* vars = 0;
  // vars = model_.getVars();
  // for (int i = 0; i < model_.get(GRB_IntAttr_NumVars); ++i)
  // {
  //   model_.remove(vars[i]);
  // }

  model_.reset();  // Note that this function, only by itself, does NOT remove vars or constraints

  ////////////////////////////
  ////////////////////////////

  ////////////////////////////  ADD CONSTRAINTS
  ////////////////////////////

  for (size_t i = 0; i < pointsA.cols(); i++)
  {
    model_.addConstr(n0_ * pointsA(0, i) + n1_ * pointsA(1, i) + n2_ * pointsA(2, i) + d_ >= epsilon_);  // n'xA+d >=
                                                                                                         // epsilon
  }

  for (size_t i = 0; i < pointsB.cols(); i++)
  {
    model_.addConstr(n0_ * pointsB(0, i) + n1_ * pointsB(1, i) + n2_ * pointsB(2, i) + d_ <= -epsilon_);  // n'xB+d
                                                                                                          // <=-epsilon
  }
  ////////////////////////////
  ////////////////////////////

  ////////////////////////////  ADD OBJECTIVE
  ////////////////////////////
  model_.setObjective(weight_n0_ * n0_ + weight_n1_ * n1_ + weight_n2_ * n2_, GRB_MINIMIZE);
  ////////////////////////////
  ////////////////////////////

  model_.update();  // needed due to the lazy evaluation
  // model_.write("/home/jtorde/Desktop/ws/src/mader/model.lp");
  model_.optimize();

  int optimstatus = model_.get(GRB_IntAttr_Status);

  num_of_LPs_run_++;  // Now num_of_LPs_run_ counts also the last LP run
  double total_time_us =
      (std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start_time))
          .count();

  // https://math.stackexchange.com/a/22351
  mean_comp_time_ms_ = mean_comp_time_ms_ + (total_time_us / 1e3 - mean_comp_time_ms_) / num_of_LPs_run_;

  std::cout << "total_time_us LP =" << total_time_us << "us" << std::endl;
  std::cout << "mean comp time LP =" << mean_comp_time_ms_ * 1000 << "us" << std::endl;
  std::cout << "num_of_LPs_run_ LP =" << num_of_LPs_run_ << std::endl;

  // int number_of_stored_solutions = model_.get(GRB_IntAttr_SolCount);
  // || optimstatus == GRB_TIME_LIMIT ||
  //       optimstatus == GRB_USER_OBJ_LIMIT ||                                    ///////////////
  //       optimstatus == GRB_ITERATION_LIMIT || optimstatus == GRB_NODE_LIMIT ||  ///////////////
  //       optimstatus == GRB_SOLUTION_LIMIT) &&
  //      number_of_stored_solutions > 0
  if (optimstatus == GRB_OPTIMAL)
  {
    solutionN(0) = n0_.get(GRB_DoubleAttr_X);
    solutionN(1) = n1_.get(GRB_DoubleAttr_X);
    solutionN(2) = n2_.get(GRB_DoubleAttr_X);
    solutionD = d_.get(GRB_DoubleAttr_X);
    return true;
  }
  else
  {
    std::cout << "Gurobi (LP) failed to find a solution" << std::endl;
    return false;
  }
};

}  // namespace separator
