cmake_minimum_required(VERSION 3.17)
project(Temperature)

set(CMAKE_CXX_STANDARD 20)

# Найдём VTK
find_package(VTK REQUIRED)
# ... и подключим её со всеми её компонентами её нативными средствами
include(${VTK_USE_FILE})

add_executable(Temperature main.cpp Matrices/CSRMatrix.h Matrices/Triplet.h Operators/partialTimeDerivative.h Field/rectangularField2D.h Operators/Laplace.h Solvers/Gauss-Seidel.h Overloads/Operators.h Overloads/Functions.h BoundaryConditions/DirihletBC.h BoundaryConditions/DirihletBC.cpp Operators/Laplace.cpp Operators/partialTimeDerivative.cpp Field/rectangularField2D.cpp BoundaryConditions/NeumannBC.cpp BoundaryConditions/NeumannBC.h BoundaryConditions/BC.h)
target_link_libraries(Temperature ${VTK_LIBRARIES})