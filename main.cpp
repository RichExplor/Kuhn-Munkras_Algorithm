/**
 * @file main.cpp
 * @author Guo Feng (1637850405@qq.com)
 * @brief 
 * @version 0.1
 * @date 2024-06-23
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#include "kuhn_munkras_algorithm.hpp"
#include <thread>
#include <memory>

int main(int agrc, char** argv)
{
    Eigen::MatrixXd cost_matrix = Eigen::MatrixXd::Zero(2, 3);
    cost_matrix << 0.0, 0.1, std::numeric_limits<double>::infinity(), 0.1, 0.1, 2.1;
    std::vector<std::pair<double, std::pair<int, int>>> match_results;
    slam::KuhnMunkrasAlgorithm KM_Algorithm(cost_matrix, false);
    KM_Algorithm.FindPerfectMatchResult(match_results);
    for (const auto& match_result : match_results)
    {
        std::cout << "cost: " << match_result.first << ", index: " << match_result.second.first << ", index: " << match_result.second.second << std::endl;
    }
    return 0;
}
