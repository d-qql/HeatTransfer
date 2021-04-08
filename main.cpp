#include <iostream>

#include "Field/rectangularField2D.h"

double start(double x, double y){

    return 600;
}
double alpha(double x, double y, double t){
    return -48.;
}
double beta(double x, double y, double t){
    return 0.01;
}
double f(double x, double y, double t){
    return 0.01*300;
}


int main() {
    RectField2D F = RectField2D(50, 50, 2);
    F.setStarterValues(start);
    F.doTimeSteps(alpha, beta, f, 5, 10, 1500, 1e-3);

    return 0;
}
