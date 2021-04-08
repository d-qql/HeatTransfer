//
// Created by d-qql on 06.04.2021.
//

#ifndef TEMPERATURE_NEUMANNBC_H
#define TEMPERATURE_NEUMANNBC_H
#include <tuple>
#include <vector>
#include <set>
#include <functional>
#include "../Matrices/Triplet.h"

class RectField2D;
std::set<Triplet> NeumannBC(const RectField2D& F, double t, const std::function<double (double x, double y, double t)>& alpha);


#endif //TEMPERATURE_NEUMANNBC_H
