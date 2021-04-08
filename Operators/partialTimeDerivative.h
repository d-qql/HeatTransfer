//
// Created by d-qql on 05.03.2021.
//

#ifndef TEMPERATURE_PARTIALTIMEDERIVATIVE_H
#define TEMPERATURE_PARTIALTIMEDERIVATIVE_H

#include <tuple>
#include <vector>
#include <set>
#include "../Matrices/Triplet.h"


class RectField2D;
std::tuple<std::set<Triplet>, std::vector<double>> pDt(const RectField2D& F, double dt);

#endif //TEMPERATURE_PARTIALTIMEDERIVATIVE_H
