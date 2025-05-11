/**
 * @file kuhn_munkras_algorithm.cpp
 * @author Guo Feng (1637850405@qq.com)
 * @brief 
 * @version 0.1
 * @date 2024-06-16
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include "kuhn_munkras_algorithm.hpp"

namespace slam
{

constexpr double MAX_MATCH_THRES = 1e8;
constexpr double MIN_MATCH_THRES = 1e-8;

KuhnMunkrasAlgorithm::KuhnMunkrasAlgorithm(const double& gate, const bool& is_max_cost_mode) {
  // 设置最大或者最小损失模式
  max_cost_mode_ = is_max_cost_mode ? 1. : -1.;
  gate_          = gate;
  is_transpose_  = false;
}

void KuhnMunkrasAlgorithm::SetCostMatrix(const Eigen::MatrixXd& cost_matrix) {
  Eigen::MatrixXd cost_matrix_bk;
  if (cost_matrix.rows() > cost_matrix.cols()) {
    cost_matrix_bk = cost_matrix.transpose();
    is_transpose_  = true;
  } else {
    cost_matrix_bk = cost_matrix;
    is_transpose_  = false;
  }

  // 根据损失权重矩阵初始化内部权重矩阵，考虑元素为0的情况
  row_                     = cost_matrix_bk.rows();
  col_                     = cost_matrix_bk.cols();
  const int max_dimensions = std::max(row_, col_);

  cost_matrix_ = Eigen::MatrixXd::Zero(max_dimensions, max_dimensions);
  for (int row = 0; row < row_; ++row) {
    for (int col = 0; col < col_; ++col) {
      cost_matrix_(row, col) = max_cost_mode_ * cost_matrix_bk(row, col);
      if (std::abs(cost_matrix_bk(row, col)) < MIN_MATCH_THRES) {
        cost_matrix_(row, col) = max_cost_mode_ * MIN_MATCH_THRES;
      }
      if (std::abs(cost_matrix_bk(row, col)) == INFINITY) {
        cost_matrix_(row, col) = max_cost_mode_ * MAX_MATCH_THRES;
      }
    }
  }

  // 左-右匹配对权重与index映射（可以根据取最大或者最小确定是损失最小或最大）
  left_match_index_  = std::vector<std::pair<double, int>>(max_dimensions, std::pair<double, int>(0.0, -1));
  right_match_index_ = std::vector<std::pair<double, int>>(max_dimensions, std::pair<double, int>(0.0, -1));
  {
    uint32_t index = 0;
    std::for_each(left_match_index_.begin(), left_match_index_.end(), [this, &index](std::pair<double, int>& elem) {
      Eigen::MatrixXd cost_matrix_row = cost_matrix_.row(index++);
      elem.first                      = cost_matrix_row.maxCoeff();
    });
  }

  // 辅助标志，用于标杆更新
  s_.resize(max_dimensions, false);
  t_.resize(max_dimensions, false);
  slack_.resize(max_dimensions, MAX_MATCH_THRES);
}

void KuhnMunkrasAlgorithm::KuhnMunkrasDfs() {
  // 基于每行最小/最大损失遍历index
  for (uint32_t left_index = 0; left_index < left_match_index_.size(); ++left_index) {
    const auto& max_cost = left_match_index_[left_index].first;
    if (std::abs(max_cost) < MIN_MATCH_THRES)
      continue;

    while (true) {
      std::generate(slack_.begin(), slack_.end(), []() { return MAX_MATCH_THRES; });
      std::generate(s_.begin(), s_.end(), []() { return false; });
      std::generate(t_.begin(), t_.end(), []() { return false; });
      // 寻找增广路径
      if (true == FindPath(left_index)) {
        break;
      } else {
        // 更新标杆
        double inner_delta = MAX_MATCH_THRES;
        for (uint32_t i = 0; i < t_.size(); ++i) {
          if (false == t_[i])
            inner_delta = std::min(inner_delta, slack_[i]);
        }

        for (uint32_t i = 0; i < t_.size(); ++i) {
          if (true == s_[i])
            left_match_index_[i].first -= inner_delta;
          if (true == t_[i])
            right_match_index_[i].first += inner_delta;
        }
      }
    }
  }
}

bool KuhnMunkrasAlgorithm::FindPath(const int& left_index) {
  bool result                  = false;
  s_[left_index]               = true;
  const auto& left_cost_matrix = cost_matrix_.row(left_index);

  for (uint32_t right_index = 0; right_index < right_match_index_.size(); ++right_index) {
    double all_cost =
        left_match_index_[left_index].first + right_match_index_[right_index].first - left_cost_matrix[right_index];
    if (false == t_[right_index] && std::abs(all_cost) < MIN_MATCH_THRES) {
      t_[right_index] = true;
      if (right_match_index_[right_index].second < 0 || FindPath(right_match_index_[right_index].second)) {
        right_match_index_[right_index].second = left_index;
        result                                 = true;
        break;
      }
    } else {
      slack_[right_index] = std::min(slack_[right_index], all_cost);
    }
  }
  return result;
}

bool KuhnMunkrasAlgorithm::CostMatrixCheck(const Eigen::MatrixXd& cost_matrix) {
  if (row_ == 0 || col_ == 0)
    return false;

  for (int i = 0; i < row_; ++i) {
    for (int j = i + 1; j < row_; ++j) {
      if (cost_matrix.row(i) == cost_matrix.row(j))
        return false;
    }
  }

  for (int i = 0; i < col_; ++i) {
    for (int j = i + 1; j < col_; ++j) {
      if (cost_matrix.col(i) == cost_matrix.col(j))
        return false;
    }
  }

  return true;
}

bool KuhnMunkrasAlgorithm::FindPerfectMatchResult(AssociateResult& match_results) {
  // 损失矩阵校验，如果存在相同行或者列，不予匹配
  if (false == CostMatrixCheck(cost_matrix_)) {
    return false;
  }

  // 基于深搜的KM匹配
  KuhnMunkrasDfs();

  // 构建最大匹配对,剔除不符合实际情况的匹配对
  std::vector<bool> used_src(is_transpose_ ? col_ : row_, false);
  std::vector<bool> used_tar(is_transpose_ ? row_ : col_, false);
  for (size_t tar_index = 0; tar_index < right_match_index_.size(); ++tar_index) {
    int src_index = right_match_index_.at(tar_index).second;
    if (src_index < 0 || std::abs(cost_matrix_(src_index, tar_index)) < MIN_MATCH_THRES) {
      continue;
    }

    double cost        = cost_matrix_(src_index, tar_index);
    bool   match_valid = max_cost_mode_ > 0 ? std::abs(cost) > gate_ : std::abs(cost) < gate_;
    if (match_valid) {
      if (!is_transpose_) {
        match_results.assignment.emplace_back(std::make_tuple(max_cost_mode_ * cost, src_index, tar_index));
        used_src.at(src_index) = true;
        used_tar.at(tar_index) = true;
      } else {
        match_results.assignment.emplace_back(std::make_tuple(max_cost_mode_ * cost, tar_index, src_index));
        used_src.at(tar_index) = true;
        used_tar.at(src_index) = true;
      }
    }
  }

  // 构建未关联的meas与track
  for (size_t src_idx = 0; src_idx < used_src.size(); ++src_idx) {
    if (!used_src.at(src_idx)) {
      match_results.unassignment_meas.emplace_back(src_idx);
    }
  }

  for (size_t tar_idx = 0; tar_idx < used_tar.size(); ++tar_idx) {
    if (!used_tar.at(tar_idx)) {
      match_results.unassignment_track.emplace_back(tar_idx);
    }
  }

  return !match_results.assignment.empty();
}

}  // namespace slam