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
#include <iostream>
#include <Eigen/Geometry>

namespace slam
{

class KuhnMunkrasAlgorithm
{
public:
    KuhnMunkrasAlgorithm(const Eigen::MatrixXd& cost_matrix, const bool& is_max_cost_mode);
    ~KuhnMunkrasAlgorithm() = default;

    bool FindPerfectMatchResult(std::vector<std::pair<double, std::pair<int, int>>>& match_results);

private:
    void KuhnMunkrasDfs();
    bool FindPath(const int& left_index);
    bool CostMatrixCheck(const Eigen::MatrixXd& cost_matrix);

private:
    Eigen::MatrixXd cost_matrix_;
    std::vector<std::pair<double, int>> left_match_index_;
    std::vector<std::pair<double, int>> right_match_index_;

    std::vector<double> slack_;
    std::vector<bool> s_;
    std::vector<bool> t_;

    double max_cost_mode_;
};

}  // namespace slam