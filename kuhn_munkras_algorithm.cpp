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

KuhnMunkrasAlgorithm::KuhnMunkrasAlgorithm(const Eigen::MatrixXd& cost_matrix)
{
    int max_dimensions = std::max(cost_matrix.cols, cost_matrix.rows);
    cost_matrix_ = Eigen::Matrix<double, max_dimensions, max_dimensions>::Zero();
    for (int row = 0; row < cost_matrix.rows; ++row)
    {
        for (int col = 0; col < cost_matrix.cols; ++col)
        {
            cost_matrix_[row][col] = cost_matrix[row][col];
        }
    }

    left_match_index_ = std::vector<std::pair<double, int>>(max_dimensions, std::pair<double, int>(0.0, -1));
    right_match_index_ = std::vector<std::pair<double, int>>(max_dimensions, std::pair<double, int>(0.0, -1));

    {
        uint32_t index = 0;
        std::for_each(left_match_index_.begin(), left_match_index_.end(),
            [this, &index])(std::pair<double, int>& elem) {
                double max_cost = cost_matrix_.at[index++].maxCoeff();
                elem.first = max_cost;
            }
    }
}

bool KuhnMunkrasAlgorithm::FindPerfectMatchResult()
{
    return true;
}

void KuhnMunkrasAlgorithm::KuhnMunkrasDfs()
{

}

bool KuhnMunkrasAlgorithm::FindPath(const int& index)
{

}

}  // namespace slam