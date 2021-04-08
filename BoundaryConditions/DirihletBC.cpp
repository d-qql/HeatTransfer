//
// Created by d-qql on 12.03.2021.
//
#include "DirihletBC.h"
#include "../Field/rectangularField2D.h"

std::set<Triplet> DirihletBC(const RectField2D& F, double t, const std::function<double (double x, double y, double t)>& beta){

    size_t H = F.getH();
    size_t W = F.getW();
    double stepX = F.getStepX();
    double stepY = F.getStepY();
    std::set<Triplet> Mpart; // matrix part of BC
    //std::vector<double> Bpart(H*W); //vector part of BC

    //i==0 || i == H-1
    for( size_t j = 0; j < W; ++j) {
    //    Bpart[j] = T;
        Mpart.insert({j, j, beta(j * stepX, 0, t)});
    //    Bpart[(H - 1) * W + j] = T;
        Mpart.insert({(H - 1) * W + j, (H - 1) * W +  j, beta(j * stepX, static_cast<double>(H - 1) * stepY, t)});
    }
    //j==0 || j == W - 1
    for( size_t i = 1; i < H - 1; ++i) {
    //    Bpart[i * W] = T;
        Mpart.insert({i * W, i * W, beta(0, i * stepY, t)});
    //    Bpart[i * W + (W - 1)] = T;
        Mpart.insert({i * W + (W - 1), i * W + (W - 1), beta(static_cast<double>(W - 1) * stepX, i * stepY, t)});
    }
    return Mpart;
}