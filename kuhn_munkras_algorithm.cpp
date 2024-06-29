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

const double& MAX_MATCH_THRES = 1E6;
const double& MIN_MATCH_THRES = 1E-6;

KuhnMunkrasAlgorithm::KuhnMunkrasAlgorithm(const Eigen::MatrixXd& cost_matrix, const bool& is_max_cost_mode)
{
    // 设置最大或者最小损失模式
    max_cost_mode_ = is_max_cost_mode ? 1. : -1.;

    // 根据损失权重矩阵初始化内部权重矩阵，考虑元素为0的情况
    const int rows = cost_matrix.rows();
    const int cols = cost_matrix.cols();
    const int max_dimensions = std::max(rows, cols);

    cost_matrix_ = Eigen::MatrixXd::Zero(max_dimensions, max_dimensions);
    // cost_matrix_.block(0, 0, rows, cols) = max_cost_mode_ * cost_matrix;
    for (int row = 0; row < rows; ++row)
    {
        for (int col = 0; col < cols; ++col)
        {
            cost_matrix_(row, col) = max_cost_mode_ * cost_matrix(row, col);
            if (std::fabs(cost_matrix(row, col)) < MIN_MATCH_THRES)
            {
                cost_matrix_(row, col) = max_cost_mode_ * MIN_MATCH_THRES;
            }
            if (std::fabs(cost_matrix(row, col)) == INFINITY)
            {
                cost_matrix_(row, col) = max_cost_mode_ * MAX_MATCH_THRES;
            }
        }
    }

    // 左-右匹配对权重与index映射（可以根据取最大或者最小确定是损失最小或最大）
    left_match_index_ = std::vector<std::pair<double, int>>(max_dimensions, std::pair<double, int>(0.0, -1));
    right_match_index_ = std::vector<std::pair<double, int>>(max_dimensions, std::pair<double, int>(0.0, -1));
    {
        uint32_t index = 0;
        std::for_each(left_match_index_.begin(), left_match_index_.end(),
            [this, &index](std::pair<double, int>& elem) {
                Eigen::MatrixXd cost_matrix_row = cost_matrix_.row(index++);
                elem.first = cost_matrix_row.maxCoeff();
            });
    }

    // 辅助标志，用于标杆更新
    s_.resize(max_dimensions, false);
    t_.resize(max_dimensions, false);
    slack_.resize(max_dimensions, MAX_MATCH_THRES);
}

void KuhnMunkrasAlgorithm::KuhnMunkrasDfs()
{
    // 基于每行最小/最大损失遍历index
    for (uint32_t left_index = 0; left_index < left_match_index_.size(); ++left_index)
    {
        const auto& max_cost = left_match_index_[left_index].first;
        if (std::fabs(max_cost) < MIN_MATCH_THRES)  continue;

        while (true)
        {
            std::generate(slack_.begin(), slack_.end(), []() { return MAX_MATCH_THRES; });
            std::generate(s_.begin(), s_.end(), []() { return false; });
            std::generate(t_.begin(), t_.end(), []() { return false; });
            // 寻找增广路径
            if (true == FindPath(left_index))
            {
                break;
            }
            else
            {
                // 更新标杆
                double inner_delta = MAX_MATCH_THRES;
                for (uint32_t i = 0; i < t_.size(); ++i)
                {
                    if (false == t_[i]) inner_delta = std::min(inner_delta, slack_[i]);
                }

                for (uint32_t i = 0; i < t_.size(); ++i)
                {
                    if (true == s_[i]) left_match_index_[i].first -= inner_delta;
                    if (true == t_[i]) right_match_index_[i].first += inner_delta;
                }
            }
        }
    }
}

bool KuhnMunkrasAlgorithm::FindPath(const int& left_index)
{
    bool result = false;
    s_[left_index] = true;
    const auto& left_cost_matrix = cost_matrix_.row(left_index);

    for (uint32_t right_index = 0; right_index < right_match_index_.size(); ++right_index)
    {
        double all_cost = left_match_index_[left_index].first + right_match_index_[right_index].first - left_cost_matrix[right_index];
        if (false == t_[right_index] && std::fabs(all_cost) < MIN_MATCH_THRES)
        {
            t_[right_index] = true;
            if (right_match_index_[right_index].second < 0 || FindPath(right_match_index_[right_index].second))
            {
                right_match_index_[right_index].second = left_index;
                result = true;
                break;
            }
        }
        else
        {
            slack_[right_index] = std::min(slack_[right_index], all_cost);
        }
    }
    return result;
}

bool KuhnMunkrasAlgorithm::CostMatrixCheck(const Eigen::MatrixXd& cost_matrix)
{
    const int rows = cost_matrix.rows();
    const int cols = cost_matrix.cols();

    if (rows == 0 || cols == 0) return false;

    for (int i = 0; i < rows; ++i)
    {
        for (int j = i + 1; j < rows; ++j)
        {
            if (cost_matrix.row(i) == cost_matrix.row(j)) return false;
        }
    }

    for (int i = 0; i < cols; ++i)
    {
        for (int j = i + 1; j < cols; ++j)
        {
            if (cost_matrix.col(i) == cost_matrix.col(j)) return false;
        }
    }

    return true;
}

bool KuhnMunkrasAlgorithm::FindPerfectMatchResult(std::vector<std::pair<double, std::pair<int, int>>>& match_results)
{
    // 损失矩阵校验，如果存在相同行或者列，不予匹配
    if (false == CostMatrixCheck(cost_matrix_))
    {
        return false;
    }

    // 基于深搜的KM匹配
    KuhnMunkrasDfs();

    // 构建最大匹配对
    for (uint32_t index = 0; index < right_match_index_.size(); ++index)
    {
        const int& src_index = right_match_index_[index].second;
        if (src_index < 0) continue;

        const double& cost = cost_matrix_(src_index, index);
        if (std::fabs(cost) < MIN_MATCH_THRES)  continue;

        match_results.emplace_back(max_cost_mode_ * cost, std::make_pair(src_index, index));
    }

    return !match_results.empty();
}

}  // namespace slam