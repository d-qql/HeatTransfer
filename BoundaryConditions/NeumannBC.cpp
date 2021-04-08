//
// Created by d-qql on 06.04.2021.
//

#include "NeumannBC.h"
#include "../Field/rectangularField2D.h"

std::set<Triplet> NeumannBC(const RectField2D& F, double t, const std::function<double (double x, double y, double t)>& alpha){

    size_t H = F.getH();
    size_t W = F.getW();
    double stepX = F.getStepX();
    double stepY = F.getStepY();

    std::set<Triplet> Mpart; // matrix part of BC
    //std::vector<double> Bpart(H*W); //vector part of BC

    //i==0 || i == H-1
    for( size_t j = 1; j < W - 1; ++j) {
        //    Bpart[j] = T;
        Mpart.insert({j, j, -alpha(j * stepX, 0, t) / stepY}); //граничные элементы
        Mpart.insert({j, W + j, alpha(j * stepX, 0, t)/stepY}); //элементы под ними
        //    Bpart[(H - 1) * W + j] = T;
        Mpart.insert({(H - 1) * W + j, (H - 1) * W +  j, -alpha(j * stepX, static_cast<double>(H - 1) * stepY, t)/stepY}); //граничные элементы
        Mpart.insert({(H - 1) * W + j, (H - 2) * W +  j, alpha(j * stepX, static_cast<double>(H - 1) * stepY, t)/stepY});    //элементы над ними
    }
    //j==0 || j == W - 1
    for( size_t i = 1; i < H - 1; ++i) {
        //    Bpart[i * W] = T;
        Mpart.insert({i * W, i * W, -alpha(0, i * stepY, t)/stepX});
        Mpart.insert({i * W, i * W + 1, alpha(0, i * stepY, t)/stepX});
        //    Bpart[i * W + (W - 1)] = T;
        Mpart.insert({i * W + (W - 1), i * W + (W - 1), -alpha(static_cast<double>(W - 1) * stepX, i * stepY, t)/stepX});
        Mpart.insert({i * W + (W - 1), i * W + (W - 2), alpha(static_cast<double>(W - 1) * stepX,i * stepY, t)/stepX});

    }

    //углы
    Mpart.insert({0, 0, 1.});
    Mpart.insert({W - 1, W - 1, 1.});
    Mpart.insert({(H - 1) * W, (H - 1) * W, 1.});
    Mpart.insert({(H - 1) * W + W - 1, (H - 1) * W + W - 1, 1.});

    return Mpart;
}