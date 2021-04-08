//
// Created by d-qql on 06.03.2021.
//

#ifndef TEMPERATURE_LAPLACE_H
#define TEMPERATURE_LAPLACE_H

#include <vector>
#include <set>
#include <cmath>
#include "../Matrices/Triplet.h"


class RectField2D;
std::set<Triplet> Laplace(const RectField2D& F);

#endif //TEMPERATURE_LAPLACE_H
