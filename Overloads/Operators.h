//
// Created by d-qql on 11.03.2021.
//

#ifndef TEMPERATURE_OPERATORS_H
#define TEMPERATURE_OPERATORS_H
#include <vector>
#include <ostream>
#include "../Matrices/CSRMatrix.h"
std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b){
    std::vector<double> res(a.size());
    for(size_t i = 0; i < a.size(); ++i) res[i] = a[i] - b[i];
    return std::move(res);
}

std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b){
    std::vector<double> res(a.size());
    for(size_t i = 0; i < a.size(); ++i) res[i] = a[i] + b[i];
    return std::move(res);
}
std::ostream& operator<<(std::ostream& os, const std::vector<double>& b){
    os<<"( ";
    for(auto i: b){
        os<<i<<" ";
    }
    os<<")\n";
    return os;
}
std::ostream& operator<<(std::ostream& os, const CSR &Matrix){
    for(size_t i = 0; i < Matrix.sizeH(); ++i){
        for(size_t j = 0; j < Matrix.sizeW(); ++j){
            os<<Matrix(i, j)<<" ";
        }
        os<<std::endl;
    }
    return os;
}

#endif //TEMPERATURE_OPERATORS_H
