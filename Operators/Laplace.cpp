//
// Created by d-qql on 12.03.2021.
//
#include "Laplace.h"
#include "../Field/rectangularField2D.h"

std::set<Triplet> Laplace(const RectField2D& F){

    size_t H = F.getH();
    size_t W = F.getW();
    double stepX = F.getStepX();
    double stepY = F.getStepY();

    std::set<Triplet> Mpart; // matrix part of derivative

    for( size_t i = 1; i < H - 1; ++i){
        for( size_t j = 1; j < W - 1; ++j) {
            Mpart.insert({i * W + j, i * W + (j - 1),  std::pow(1. / stepY, 2)});
            Mpart.insert({i * W + j, i * W + (j + 1),  std::pow(1. / stepY, 2)});

            Mpart.insert({i * W + j, i * W + j, -2. * (std::pow(1. / stepX, 2) + std::pow(1. / stepY, 2))});

            Mpart.insert({i * W + j, (i - 1) * W + j,  std::pow(1. / stepX, 2)});
            Mpart.insert({i * W + j, (i + 1) * W + j,  std::pow(1. / stepX, 2)});
        }
    }
    /*
    // i == 0 || i == H - 1:
    for( size_t j = 1; j < W - 1; ++j) {
        Mpart.insert({j, j - 1,  std::pow(1. / stepY, 2)}); // i, j - 1
        Mpart.insert({j, j, -2. * (std::pow(1. / stepX, 2) + std::pow(1. / stepY, 2))}); // i , j
        Mpart.insert({j, j + 1,  std::pow(1. / stepY, 2)}); //i, j + 1
        Mpart.insert({ j, (1 * W) + j,  std::pow(1. / stepX, 2)}); // i + 1, j

        Mpart.insert({(H - 1) * W + j, (H - 1) * W + j - 1,  std::pow(1. / stepY, 2)}); // j - 1
        Mpart.insert({(H - 1) * W + j, (H - 1) * W + j, -2. * (std::pow(1. / stepX, 2) + std::pow(1. / stepY, 2))}); // i , j
        Mpart.insert({(H - 1) * W + j, (H - 1) * W + j + 1,  std::pow(1. / stepY, 2)}); // j + 1
        Mpart.insert({(H - 1) * W + j, (H - 1 - 1) * W + j,  std::pow(1. / stepX, 2)}); // i - 1
    }

    // j == 0 || j == W - 1
    for( size_t i = 1; i < H - 1; ++i) {
        Mpart.insert({i * W, (i - 1) * W,  std::pow(1. / stepX, 2)}); // i - 1, j
        Mpart.insert({i * W, i * W, -2. * (std::pow(1. / stepX, 2) + std::pow(1. / stepY, 2))}); // i , j
        Mpart.insert({i * W, (i + 1) * W,  std::pow(1. / stepX, 2)}); //i + 1, j
        Mpart.insert({i * W, i * W + 1,  std::pow(1. / stepY, 2)}); // j + 1

        Mpart.insert({i * W + (W - 1), (i - 1) * W + (W - 1),  std::pow(1. / stepX, 2)}); // i - 1, j
        Mpart.insert({i * W + (W - 1), i * W + (W - 1), -2. * (std::pow(1. / stepX, 2) + std::pow(1. / stepY, 2))}); // i , j
        Mpart.insert({i * W + (W - 1), (i + 1) * W + (W - 1),  std::pow(1. / stepX, 2)}); // i + 1, j
        Mpart.insert({i * W + (W - 1), i * W + (W - 2),  std::pow(1. / stepY, 2)}); // j - 1
    }

    //i == 0 && j == 0:
    Mpart.insert({0, 0, -2. * (std::pow(1. / stepX, 2) + std::pow(1. / stepY, 2))});
    Mpart.insert({0, 1, std::pow(1. / stepY, 2)}); //j+1
    Mpart.insert({0, W, std::pow(1. / stepX, 2)}); //i+1

    //i == 0 && j == W - 1:
    Mpart.insert({W - 1, (W - 1), -2. * (std::pow(1. / stepX, 2) + std::pow(1. / stepY, 2))});
    Mpart.insert({W - 1, (W - 1 - 1), std::pow(1. / stepY, 2)}); //j-1
    Mpart.insert({W - 1, W + (W - 1), std::pow(1. / stepX, 2)}); //i+1

    //i == H - 1 && j == 0:
    Mpart.insert({(H - 1) * W, (H - 1) * W + 0, -2. * (std::pow(1. / stepX, 2) + std::pow(1. / stepY, 2))});
    Mpart.insert({(H - 1) * W, (H - 1) * W + 1, std::pow(1. / stepY, 2)}); //j+1
    Mpart.insert({(H - 1) * W, (H - 1 - 1) * W + 0, std::pow(1. / stepX, 2)}); //i-1

    //i == H - 1 && j == W - 1:
    Mpart.insert({(H - 1) * W + (W - 1), (H - 1) * W + W - 1, -2. * (std::pow(1. / stepX, 2) + std::pow(1. / stepY, 2))});
    Mpart.insert({(H - 1) * W + (W - 1), (H - 1) * W + W - 2, std::pow(1. / stepY, 2)}); //j-1
    Mpart.insert({(H - 1) * W + (W - 1), (H - 1 - 1) * W + W - 1, std::pow(1. / stepY, 2)}); //i-1
*/

    return Mpart;
}