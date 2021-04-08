//
// Created by d-qql on 11.03.2021.
//

#ifndef TEMPERATURE_FUNCTIONS_H
#define TEMPERATURE_FUNCTIONS_H
#include <vector>
#include <cmath>

double norm(const std::vector<double>& vec){
    double sum = 0;
    for(const auto& elm : vec){
        sum += std::pow(elm, 2);
    }
    return std::sqrt(sum);
}

#endif //TEMPERATURE_FUNCTIONS_H
