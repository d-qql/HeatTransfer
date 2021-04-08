//
// Created by d-qql on 11.03.2021.
//

#ifndef TEMPERATURE_DIRIHLETBC_H
#define TEMPERATURE_DIRIHLETBC_H

#include <tuple>
#include <vector>
#include <set>
#include <functional>
#include "../Matrices/Triplet.h"

class RectField2D;
std::set<Triplet> DirihletBC(const RectField2D& F, double t, const std::function<double (double x, double y, double t)>& beta);

#endif //TEMPERATURE_DIRIHLETBC_H
