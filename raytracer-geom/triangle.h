#pragma once

#include "vector.h"

#include <cstddef>
#include <stdexcept>
class Triangle {
public:
    Triangle(const Vector& a, const Vector& b, const Vector& c) : a_(a), b_(b), c_(c) {
    }

    const Vector& operator[](size_t ind) const {
        if (ind == 0) {
            return a_;
        }
        if (ind == 1) {
            return b_;
        }
        if (ind == 2) {
            return c_;
        }
        throw std::runtime_error("Index should be 0/1/2");
    }
    double Area() const {
        return Length(CrossProduct((b_ - a_), (c_ - a_))) * 1 / 2;
    }

private:
    Vector a_;
    Vector b_;
    Vector c_;
};
