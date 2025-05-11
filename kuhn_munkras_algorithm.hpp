/**
 * @file kuhn_munkras_algorithm.hpp
 * @author Guo Feng (1637850405@qq.com)
 * @brief 
 * @version 0.1
 * @date 2024-06-16
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#pragma once

#include <set>
#include <vector>
#include <tuple>
#include <iostream>
#include <Eigen/Geometry>

namespace slam
{

struct AssociateResult {
  std::vector<std::tuple<double, int, int>> assignemnt;
  std::vector<int>                          unassignment_meas;
  std::vector<int>                          unassignment_track;
};

class KuhnMunkrasAlgorithm {
public:
  KuhnMunkrasAlgorithm(const double& gate = 3.0, const bool& is_max_cost_mode = false);
  ~KuhnMunkrasAlgorithm() = default;

  void SetCostMatrix(const Eigen::MatrixXd& cost_matrix);
  bool FindPerfectMatchResult(AssociateResult& match_results);

private:
  void KuhnMunkrasDfs();
  bool FindPath(const int& left_index);
  bool CostMatrixCheck(const Eigen::MatrixXd& cost_matrix);

private:
  Eigen::MatrixXd                     cost_matrix_;
  std::vector<std::pair<double, int>> left_match_index_;
  std::vector<std::pair<double, int>> right_match_index_;

  std::vector<double> slack_;
  std::vector<bool>   s_;
  std::vector<bool>   t_;
  int                 row_;
  int                 col_;

  double max_cost_mode_;
  double gate_;
  bool   is_transpose_;
};

}  // namespace slam