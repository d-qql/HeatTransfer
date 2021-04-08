//
// Created by d-qql on 05.03.2021.
//

#ifndef TEMPERATURE_RECTANGULARFIELD2D_H
#define TEMPERATURE_RECTANGULARFIELD2D_H

#define eps 1e-10
#include <vector>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>

class RectField2D{
public:

    using idx_t = std::size_t;
    using val_t = double;

protected:

    idx_t H, W; //число узлов по вертикали и горизонтали
    double stepX, stepY; //шаг по пространству между узлами
    std::vector<val_t> nodes;   //вектор значений в узлах

public:

    RectField2D(const idx_t& H, const idx_t& W, const double& stepX, const double& stepY);

    RectField2D(const idx_t& H, const idx_t& W, const double& step);

    void setStarterValues(idx_t number, val_t value);

    void setStarterValues(idx_t i, idx_t j, val_t value);

    void setStarterValues(std::function<double (double, double )>);

    [[nodiscard]] std::size_t getH() const;

    [[nodiscard]] std::size_t getW() const;

    [[nodiscard]] std::size_t getLength() const;

    [[nodiscard]] double getStepX() const;

    [[nodiscard]] double getStepY() const;


    val_t operator()(const idx_t& i, const idx_t& j) const;
    val_t operator()(const idx_t& number) const;

    void doTimeSteps(const std::function<double (double x, double y, double t)>& alpha,
                     const std::function<double (double x, double y, double t)>& beta,
                     const std::function<double (double x, double y, double t)>& f,
                     double step, double a, unsigned int N, double tolerance = eps);

    // Метод отвечает за запись текущего состояния сетки в снапшот в формате VTK
    void snapshot(unsigned int snap_number);


};


#endif //TEMPERATURE_RECTANGULARFIELD2D_H
