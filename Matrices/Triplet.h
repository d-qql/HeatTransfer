//
// Created by d-qql on 05.03.2021.
//

#ifndef TEMPERATURE_TRIPLET_H
#define TEMPERATURE_TRIPLET_H


struct Triplet{

    using idx_t = std::size_t;
    using elm_t = double;

    idx_t i;
    idx_t j;
    elm_t value;

    bool operator<(Triplet const &rgh) const{
        return this->i<rgh.i || (this->i == rgh.i && this->j < rgh.j);
    }
};

#endif //TEMPERATURE_TRIPLET_H
