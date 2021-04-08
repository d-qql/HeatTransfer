//
// Created by d-qql on 05.03.2021.
//

#ifndef TEMPERATURE_CSRMATRIX_H
#define TEMPERATURE_CSRMATRIX_H

#define eps 1e-10

#include <utility>
#include <vector>
#include <set>
#include "Triplet.h"


class CSR{
public:

    using idx_t = std::size_t;
    using elm_t = double;

private:

    idx_t H, W;
    std::vector<elm_t> values;
    std::vector<idx_t> cols;
    std::vector<idx_t> rows;

    friend std::vector<elm_t> GaussSeidel(const CSR& A, const std::vector<elm_t>& b, std::vector<double> x, double tolerance);

public:
    CSR(const idx_t &H, const idx_t &W, const std::set<Triplet> &in): H(H), W(W){
        values.resize(in.size());
        cols.resize(in.size());
        rows.resize(H+1);
        int countInRow = 0;
        int currRow = 0;
        auto it = in.begin();                                   //получаем итератор на первый ненулевой элемент из set
        for(idx_t k = 0; k < in.size(); ++k){                   //проход по всем ненулевым элементам
            while(currRow < it->i){                             //Пока номера текущих обрабатываемых строк меньше номера строки текущего элемента, подлежащего вставке:
                rows[currRow + 1] = rows[currRow] + countInRow; //вставляем суммарное количество ненулевых элементов во всех предыдущих строках + в текущей
                ++currRow;                                       //переходим к заполнению следующей строки
                countInRow = 0;                                   //обнуляем счетчик вставленных в текущую строку элементов
            }                                                    //если номер текущей строки совпал с номером строки элемента, подлежащего вставке, то:
            values[k] = it->value;              //вставляем значение элемента
            cols[k] = it->j;                    //вставлем номер столбца элемента
            ++countInRow;                       //увеличиваем счетчик вставленных в текущую строку элементов
            it = next(it);                      //переходим к следующему элементу
        }

        for( ++currRow; currRow <= H; ++currRow ) rows[currRow] = in.size();  //если в конце currRow < H => имеются нулевые строки в конце матрицы
        // (необходимо заполнить их индексацию числом ненулевых элементов всей матрицы)
    }

    CSR(idx_t H, idx_t W, std::vector<elm_t> vals, std::vector<idx_t> cols, std::vector<idx_t> rows)
        : H(H), W(W), values(std::move(vals)), cols(std::move(cols)), rows(std::move(rows)){}

    [[nodiscard]] idx_t sizeH() const {
        return H;
    }
    [[nodiscard]] idx_t sizeW() const {
        return W;
    }

    elm_t operator()(idx_t const i, idx_t const j) const{
        idx_t skip = this->rows[i];
        idx_t count = this->rows[i+1] - this->rows[i];
        for (idx_t k = skip; k < skip+count; ++k) {
            if(this->cols[k] == j) return this->values[k];
        }
        return static_cast<elm_t>(0);
    }

    std::vector<elm_t> operator*(const std::vector<double>& x) const{
        std::vector<elm_t> res(this->H);
        for(idx_t i = 0; i < this->H; ++i){
            res[i] = 0.;
            for(idx_t j = this->rows[i]; j < this->rows[i+1]; ++j) res[i] += this->values[j] * x[this->cols[j]];
        }
        return res;
    }

    CSR operator*=(double k){
        for(double & value : values){
            value *= k;
        }
        return *this;
    }

    CSR operator+(const CSR &B) const {
        std::vector<idx_t> IX(this->W, -1);
        std::vector<idx_t> t_JC; //временный массив индексации столбцов для маски
        std::vector<idx_t> JC;  //массив индексации столбцов результирующей матрицы
        std::vector<idx_t> R(this->H + 1);
        std::vector<elm_t> X; //временный массив значений сложенных строк
        std::vector<elm_t> V;
        R[0] = 0;

        // слияние индексов столбцов ненулевых элементов в один массив строки i
        for(size_t i = 0; i < this->H; ++i) {
            X.resize(this->W);
            for(idx_t j = this->rows[i]; j < this->rows[i + 1]; ++j){
                if(IX[this->cols[j]] != i){
                    IX[this->cols[j]] = i;          //ставим метку, что в текущей строке столбец cols[j] учтен
                    t_JC.emplace_back(this->cols[j]);   //во временный масив столбцов записываем номер столбца
                }
                X[this->cols[j]] = this->values[j];  //записать значения строки из массива А в суммирующую строку(плотное представление строки)
            }
            for(idx_t j = B.rows[i]; j < B.rows[i + 1]; ++j){
                if(IX[B.cols[j]] != i){
                    IX[B.cols[j]] = i;
                    t_JC.emplace_back(B.cols[j]);
                }
                X[B.cols[j]] += B.values[j]; //прибавление значений из строки массива B к суммирующей строке
            }
            for(const auto& j: t_JC){
                if( std::abs(X[j]) > eps ){
                    JC.emplace_back(j);
                    V.emplace_back(X[j]);
                } else continue;
            }
            R[i + 1] = JC.size();
            t_JC.clear();
            X.clear();
        }
        CSR C(this->H, this->W, V, JC, R);
        return std::move(C);
    }

};


#endif //TEMPERATURE_CSRMATRIX_H
