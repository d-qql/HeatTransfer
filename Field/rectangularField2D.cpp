//
// Created by d-qql on 12.03.2021.
//
#include "../Overloads/Operators.h"
#include "../Solvers/Gauss-Seidel.h"
#include "rectangularField2D.h"
#include "../Operators/Laplace.h"
#include "../Operators/partialTimeDerivative.h"
#include "../BoundaryConditions/DirihletBC.h"
#include "../BoundaryConditions/NeumannBC.h"

RectField2D::RectField2D(const idx_t& H, const idx_t& W, const double& stepX, const double& stepY): H(H), W(W), stepX(stepX), stepY(stepY){
    nodes.resize(H*W);
}

RectField2D::RectField2D(const idx_t& H, const idx_t& W, const double& step): H(H), W(W), stepX(step), stepY(step){
    nodes.resize(H*W);
}

void RectField2D::setStarterValues(idx_t number, val_t value){
    nodes[number] = value;
}

void RectField2D::setStarterValues(idx_t i, idx_t j, val_t value){
    nodes[i * W + j] = value;
}

void RectField2D::setStarterValues(std::function<double (double, double )> func){
    for(size_t i = 1; i < H - 1; ++i){
        for(size_t j = 1; j < W - 1; ++j){
            nodes[i * W + j] = func(j * stepX, i * stepY);
        }
    }

}

[[nodiscard]] size_t RectField2D::getH() const{
    return this->H;
}

[[nodiscard]] size_t RectField2D::getW() const{
    return this->H;
}

[[nodiscard]] size_t RectField2D::getLength() const{
    return nodes.size();
}

[[nodiscard]] double RectField2D::getStepX() const{
    return this->stepX;
}

[[nodiscard]] double RectField2D::getStepY() const{
    return this->stepY;
}


RectField2D::val_t RectField2D::operator()(const idx_t& i, const idx_t& j) const {
    return nodes[i*W+j];
}
RectField2D::val_t RectField2D::operator()(const idx_t& number) const {
    return nodes[number];
}

void RectField2D::doTimeSteps(const std::function<double (double x, double y, double t)>& alpha,
                              const std::function<double (double x, double y, double t)>& beta,
                              const std::function<double (double x, double y, double t)>& rightPart,
                              double step, double a, unsigned int N, double tolerance){
    size_t sz;
    sz = W * H;
    std::vector<double> f(H * W);
    double t = 0; // current time
    for(unsigned int i = 0; i < N; ++i) {
        CSR A(sz, sz, Laplace(*this));
        auto[A_, b_] = pDt(*this, step);
        A *= -a;
        A = A + CSR(sz, sz, A_) + CSR(sz, sz, NeumannBC(*this, t, alpha)) + CSR(sz, sz, DirihletBC(*this, t, beta));

        for(size_t x = 0; x < W; ++x){
            b_[x] += rightPart(x * stepX, 0, t);
            b_[(H-1) * W + x] += rightPart(x * stepX, static_cast<double>(H-1) * stepY, t);
        }
        for(size_t y = 1; y < H - 1; ++y){
            b_[y * W] += rightPart(0, y * stepY, t);
            b_[y * W + (W - 1)] += rightPart(static_cast<double>(W - 1) * stepX, y * stepY, t);
        }

//          std::cout<<A;
//          std::cout<<i<<std::endl;
//          std::cout<<f+b_;
//        double sum = 0;
//        for(int k = 0; k < H; ++k){
//            for(int p = 0; p < W; ++p){
//                if (!(p == 0 && k == 0 || p == 0 && k == H - 1 || p == W - 1 && k == 0 || p == W - 1 && k == H - 1)) sum +=nodes[k * W + p];
//            }
//        }
        nodes = GaussSeidel(A, b_, nodes, tolerance);
        t += step;
        //std::cout<<nodes;
        std::cout<<i<<std::endl;
        snapshot(i);
    }
}

// Метод отвечает за запись текущего состояния сетки в снапшот в формате VTK
void RectField2D::snapshot(unsigned int snap_number) {
    // Сетка в терминах VTK
    vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
    // Точки сетки в терминах VTK
    vtkSmartPointer<vtkPoints> dumpPoints = vtkSmartPointer<vtkPoints>::New();

    // Скалярное поле на точках сетки
    auto T = vtkSmartPointer<vtkDoubleArray>::New();
    T->SetName("Temperature");

    // Обходим все точки нашей расчётной сетки
    unsigned int number = (unsigned int)nodes.size();
    for(unsigned int i = 0; i < H; i++) {
        for(unsigned int j = 0; j < W; j++) {
            // Вставляем новую точку в сетку VTK-снапшота
            dumpPoints->InsertNextPoint(i*stepX, j*stepY, 0);

            // значение скалярного поля
            T->InsertNextValue(nodes[i * W + j]);
        }
    }

    // Задаём размеры VTK-сетки (в точках, по трём осям)
    structuredGrid->SetDimensions(H, W, 1);
    // Грузим точки в сетку
    structuredGrid->SetPoints(dumpPoints);

    // Присоединяем векторное и скалярное поля к точкам
    structuredGrid->GetPointData()->AddArray(T);

    // Создаём снапшот в файле с заданным именем
    std::string fileName = "../VTK/Temperature2d-step-" + std::to_string(snap_number) + ".vts";
    vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
    writer->SetFileName(fileName.c_str());
    writer->SetInputData(structuredGrid);
    writer->Write();
}