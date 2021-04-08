//
// Created by d-qql on 11.03.2021.
//

#ifndef TEMPERATURE_GAUSS_SEIDEL_H
#define TEMPERATURE_GAUSS_SEIDEL_H

#include "../Matrices/CSRMatrix.h"
#include "../Overloads/Operators.h"
#include "../Overloads/Functions.h"

std::vector<double> GaussSeidel(const CSR& A, const std::vector<double>& b, std::vector<double> x, double tolerance = eps){
    using idx_t = typename CSR::idx_t;

    std::vector<double> r(b.size());

    std::vector<std::pair<int, double>> plotData;
    double sum;
    r = A * x - b;
    while (norm(r) > tolerance) {
        for (idx_t i = 0; i < A.H; ++i) {
            sum = 0.;
            for (idx_t j = A.rows[i]; j < A.rows[i + 1]; ++j) {
                if (i != A.cols[j]) {
                    sum += A.values[j] * x[A.cols[j]];
                }else{ continue; }
            }
            x[i] = (b[i] - sum)/A(i, i);
        }
        r = A * x - b;
    }
    return std::move(x);
}
#endif //TEMPERATURE_GAUSS_SEIDEL_H
