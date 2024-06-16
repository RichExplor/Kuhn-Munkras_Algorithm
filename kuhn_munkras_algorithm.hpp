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

#include <vector>
#include <iostream>
#include <Eigen/Geometry>

namespace slam
{
class KuhnMunkrasAlgorithm
{
public:
    KuhnMunkrasAlgorithm(const Eigen::MatrixXd& cost_matrix);
    ~KuhnMunkrasAlgorithm() = default;

    bool FindPerfectMatchResult();

private:
    void KuhnMunkrasDfs();
    bool FindPath(const int& index);

private:
    Eigen::MatrixXd cost_matrix_;
    std::vector<std::pair<double, int>> left_match_index_;
    std::vector<std::pair<double, int>> right_match_index_;

    std::vector<double> slack_;
    std::vector<bool> s_;
    std::vector<bool> t_;

};

}  // namespace slam