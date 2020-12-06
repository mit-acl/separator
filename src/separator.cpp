/* ----------------------------------------------------------------------------
 * Copyright 2020, Jesus Tordesillas Torres, Aerospace Controls Laboratory
 * Massachusetts Institute of Technology
 * All Rights Reserved
 * Authors: Jesus Tordesillas, et al.
 * See LICENSE file for the license information
 * -------------------------------------------------------------------------- */

#include <time.h>
#include <math.h>
#include <glpk.h> /* GNU GLPK linear/mixed integer solver */

#include <chrono>

#include "separator.hpp"
#include <iostream>

#define M_PI 3.14159265358979323846

namespace separator
{
Separator::Separator()  // double weight_n1, double weight_n2, double weight_n3
{
  glp_init_smcp(&params_);
  params_.msg_lev = 1;  // 1=no output.  GLP_MSG_ALL;

  weight_n1_ = 0.0;
  weight_n2_ = 0.0;
  weight_n3_ = 0.0;
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

  lp_ = glp_create_prob();
  glp_set_prob_name(lp_, "separator");
  glp_set_obj_dir(lp_, GLP_MAX);

  /* fill problem */
  glp_add_rows(lp_, pointsA.cols() + pointsB.cols());
  int row = 1;

  // See here why we can use an epsilon of 1.0:
  // http://www.joyofdata.de/blog/testing-linear-separability-linear-programming-r-glpk/
  // This also allows you to avoid checking if norm(n)==0 at the end (see degenerateSolution commented out below)
  double epsilon = 1.0;

  // n (the solution) will point to the pointsA

  // std::cout << "Using PointsA=" << std::endl;

  for (int i = 0; i < pointsA.cols(); i++)
  {
    // glp_set_row_name(lp, r, "p");
    glp_set_row_bnds(lp_, row, GLP_LO, epsilon, 0.0);  // n'xA+d>=epsilon
    row++;
  }
  //  std::cout << "Using PointsB=" << std::endl;
  for (int i = 0; i < pointsB.cols(); i++)
  {
    // glp_set_row_name(lp, r, "p");
    glp_set_row_bnds(lp_, row, GLP_UP, 0.0, -epsilon);  //<=0.0   n'xB+d <=-epsilon
    row++;
  }

  // glp_set_row_bnds(lp_, row, GLP_UP, 0.0, 3.0);  // n1+n2+n3<=3 //To prevent unbounded solutions
  row++;

  ///
  glp_add_cols(lp_, 4);

  // weights
  glp_set_col_name(lp_, 1, "n1");
  glp_set_col_bnds(lp_, 1, GLP_FR, 0.0, 0.0);  // Free
  glp_set_obj_coef(lp_, 1, weight_n1_);        // weight on n1

  glp_set_col_name(lp_, 2, "n2");
  glp_set_col_bnds(lp_, 2, GLP_FR, 0.0, 0.0);  // Free
  glp_set_obj_coef(lp_, 2, weight_n2_);        // weight on n2

  glp_set_col_name(lp_, 3, "n3");
  glp_set_col_bnds(lp_, 3, GLP_FR, 0.0, 0.0);  // Free
  glp_set_obj_coef(lp_, 3, weight_n3_);        // weight on n3

  glp_set_col_name(lp_, 4, "d");
  glp_set_col_bnds(lp_, 4, GLP_FR, 0.0, 0.0);  // Free
  glp_set_obj_coef(lp_, 4, 0.0);               // weight on n4

  int r = 1;
  row = 1;
  for (int i = 0; i < pointsA.cols(); i++)
  {
    ia_[r] = row, ja_[r] = 1, ar_[r] = pointsA(0, i); /* a[1,1] = 1 */
    r++;
    ia_[r] = row, ja_[r] = 2, ar_[r] = pointsA(1, i);  // a[1,2] = 2
    r++;
    ia_[r] = row, ja_[r] = 3, ar_[r] = pointsA(2, i);  // a[1,3] = 2
    r++;
    ia_[r] = row, ja_[r] = 4, ar_[r] = 1.0;  // a[1,4] = 1
    r++;
    row++;
  }

  for (int i = 0; i < pointsB.cols(); i++)
  {
    ia_[r] = row, ja_[r] = 1, ar_[r] = pointsB(0, i);  // a[1,1] = 1
    r++;
    ia_[r] = row, ja_[r] = 2, ar_[r] = pointsB(1, i);  // a[1,2] = 2
    r++;
    ia_[r] = row, ja_[r] = 3, ar_[r] = pointsB(2, i);  // a[1,3] = 2
    r++;
    ia_[r] = row, ja_[r] = 4, ar_[r] = 1.0;  // a[1,4] = 2
    r++;
    row++;
  }

  glp_load_matrix(lp_, r - 1, ia_, ja_,
                  ar_);  // need r-1 to substract from r++ in the last iteration of the previous loop
  // glp_write_lp(lp_, NULL, "/home/jtorde/Desktop/ws/src/faster/faster/my_model.txt");

  /* solve problem */
  // auto start = std::chrono::high_resolution_clock::now();
  glp_simplex(lp_, &params_);
  // auto elapsed = std::chrono::high_resolution_clock::now() - start;
  // std::cout << "Solving an LP took " << std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count() /
  // 1000.0
  //           << " ms" << std::endl;

  // std::cout << std::chrono::duration_cast<std::chrono::microseconds>(elapsed).count() / 1000.0 << std::endl;

  /* recover and display results */
  double z = glp_get_obj_val(lp_);
  solutionN(0) = glp_get_col_prim(lp_, 1);
  solutionN(1) = glp_get_col_prim(lp_, 2);
  solutionN(2) = glp_get_col_prim(lp_, 3);
  solutionD = glp_get_col_prim(lp_, 4);
  // std::cout << "solutionD= " << solutionD << std::endl;
  // printf("z = %g; n1 = %g; n2 = %g; n3 = %g\n", z, n1, n2, n3);

  int status = glp_get_status(lp_);

  // std::cout << "status= " << status << std::endl;

  // /*  GLP_OPT — solution is optimal;
  //   GLP_FEAS — solution is feasible;
  //   GLP_INFEAS — solution is infeasible;
  //   GLP_NOFEAS — problem has no feasible solution;
  //   GLP_UNBND — problem has unbounded solution;
  //   GLP_UNDEF — solution is undefined.*/

  // switch (status)
  // {
  //   case GLP_OPT:
  //     std::cout << "status = GLP_OPT" << std::endl;
  //     break;
  //   case GLP_FEAS:
  //     std::cout << "status = GLP_FEAS" << std::endl;
  //     break;
  //   case GLP_INFEAS:
  //     std::cout << "status = GLP_INFEAS" << std::endl;
  //     break;
  //   case GLP_NOFEAS:
  //     std::cout << "status = GLP_NOFEAS" << std::endl;
  //     break;
  //   case GLP_UNBND:
  //     std::cout << "status = GLP_UNBND" << std::endl;
  //     break;
  //   case GLP_UNDEF:
  //     std::cout << "status = GLP_UNDEF" << std::endl;
  //     break;
  //   default:
  //     std::cout << "This code doesn't exist!!" << std::endl;
  // }

  /*  if ((status != GLP_OPT) && (status != GLP_FEAS))
    {
      glp_write_lp(lp_, NULL, "/home/jtorde/Desktop/ws/src/faster/faster/my_model2.txt");
    }*/

  glp_delete_prob(lp_);
  glp_free_env();

  // std::cout << "solutionN.norm()=" << solutionN.norm() << std::endl;

  // bool degenerateSolution = (solutionN.norm() < 0.000001);  // solution is [0 0 0]
  num_of_LPs_run_++;  // Now num_of_LPs_run_ counts also the last LP run
  double total_time_us =
      (std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start_time))
          .count();

  // std::cout << "total_time_us LP =" << total_time_us << "us" << std::endl;
  // std::cout << "mean comp time LP =" << mean_comp_time_ms_ * 1000 << "us" << std::endl;
  // std::cout << "num_of_LPs_run_ LP =" << num_of_LPs_run_ << std::endl;

  // https://math.stackexchange.com/a/22351
  mean_comp_time_ms_ = mean_comp_time_ms_ + (total_time_us / 1e3 - mean_comp_time_ms_) / num_of_LPs_run_;

  if ((status == GLP_OPT || status == GLP_FEAS))  //&& !degenerateSolution
  {
    return true;
  }

  return false;
};

}  // namespace separator
