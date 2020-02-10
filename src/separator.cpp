#include <time.h>
#include <math.h>
#include <glpk.h> /* GNU GLPK linear/mixed integer solver */

#include "separator.hpp"
#include <iostream>

#define M_PI 3.14159265358979323846

namespace separator
{
Separator::Separator(double weight_n1, double weight_n2, double weight_n3)
{
  glp_init_smcp(&params_);
  params_.msg_lev = 1;  // 1=no output.  GLP_MSG_ALL-;

  weight_n1_ = weight_n1;
  weight_n2_ = weight_n2;
  weight_n3_ = weight_n3;
};

bool Separator::solveModel(Eigen::Vector3d& solution, const std::vector<Eigen::Vector3d>& pointsA,
                           const std::vector<Eigen::Vector3d>& pointsB)
{
  lp_ = glp_create_prob();
  glp_set_prob_name(lp_, "separator");
  glp_set_obj_dir(lp_, GLP_MAX);

  /* fill problem */
  glp_add_rows(lp_, pointsA.size() + pointsB.size());
  int row = 1;
  for (auto pointsA_i : pointsA)
  {
    // glp_set_row_name(lp, r, "p");
    glp_set_row_bnds(lp_, row, GLP_UP, 0.0, -1.0);  //<=-1
    row++;
  }
  for (auto pointsB_i : pointsB)
  {
    // glp_set_row_name(lp, r, "p");
    glp_set_row_bnds(lp_, row, GLP_LO, -1.0, 0.0);  //>=-1
    row++;
  }

  ///
  glp_add_cols(lp_, 3);

  glp_set_col_name(lp_, 1, "n1");
  glp_set_col_bnds(lp_, 1, GLP_FR, 0.0, 0.0);  // Free
  glp_set_obj_coef(lp_, 1, weight_n1_);        // weight on n1

  glp_set_col_name(lp_, 2, "n2");
  glp_set_col_bnds(lp_, 2, GLP_FR, 0.0, 0.0);  // Free
  glp_set_obj_coef(lp_, 2, weight_n2_);        // weight on n2

  glp_set_col_name(lp_, 3, "n3");
  glp_set_col_bnds(lp_, 3, GLP_FR, 0.0, 0.0);  // Free
  glp_set_obj_coef(lp_, 3, weight_n3_);        // weight on n3

  int r = 1;
  row = 1;
  for (auto pointsA_i : pointsA)
  {
    ia_[r] = row, ja_[r] = 1, ar_[r] = pointsA_i.x(); /* a[1,1] = 1 */
    r++;
    ia_[r] = row, ja_[r] = 2, ar_[r] = pointsA_i.y();  // a[1,2] = 2
    r++;
    ia_[r] = row, ja_[r] = 3, ar_[r] = pointsA_i.z();  // a[1,3] = 2
    r++;
    row++;
  }

  for (auto pointsB_i : pointsB)
  {
    ia_[r] = row, ja_[r] = 1, ar_[r] = pointsB_i.x();  // a[1,1] = 1
    r++;
    ia_[r] = row, ja_[r] = 2, ar_[r] = pointsB_i.y();  // a[1,2] = 2
    r++;
    ia_[r] = row, ja_[r] = 3, ar_[r] = pointsB_i.z();  // a[1,3] = 2
    r++;
    row++;
  }

  glp_load_matrix(lp_, r - 1, ia_, ja_, ar_);
  // glp_write_lp(lp_, NULL, "/home/jtorde/Desktop/ws/src/faster/faster/my_model.txt");

  /* solve problem */
  glp_simplex(lp_, &params_);
  /* recover and display results */
  double z = glp_get_obj_val(lp_);
  solution(0) = glp_get_col_prim(lp_, 1);
  solution(1) = glp_get_col_prim(lp_, 2);
  solution(2) = glp_get_col_prim(lp_, 3);
  // printf("z = %g; n1 = %g; n2 = %g; n3 = %g\n", z, n1, n2, n3);

  int status = glp_get_status(lp_);

  glp_delete_prob(lp_);
  glp_free_env();

  if (status == GLP_OPT || status == GLP_FEAS)
  {
    return true;
  }
  return false;
};

/*void Separator::deleteModel(){

};*/

}  // namespace separator
