//
// Created by d-qql on 12.03.2021.
//
#include "partialTimeDerivative.h"
#include "../Field/rectangularField2D.h"

std::tuple<std::set<Triplet>, std::vector<double>> pDt(const RectField2D& F, double dt){

    size_t H = F.getH();
    size_t W = F.getW();

    std::set<Triplet> Mpart; // matrix part of derivative
    std::vector<double> Bpart(H*W); //vector part of derivative

    for( size_t i = 1; i < H - 1; ++i){
        for( size_t j = 1; j < W - 1; ++j) {
            Bpart[i * H + j] = (F(i, j)) / dt;
            Mpart.insert({i * W + j, i * W + j, 1. / dt});
        }
    }
    return std::tuple<std::set<Triplet>, std::vector<double>>(Mpart, Bpart);
}