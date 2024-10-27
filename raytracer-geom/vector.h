#pragma once

#include <array>
#include <cstddef>
#include <cmath>
#include <stdexcept>

class Vector {
public:
    Vector() : data_{0, 0, 0} {
    }
    Vector(double x, double y, double z) {
        data_[0] = x;
        data_[1] = y;
        data_[2] = z;
    }

    double& operator[](size_t ind) {
        return data_[ind];
    }
    double operator[](size_t ind) const {
        return data_[ind];
    }
    double operator*(const Vector& other) const {
        return data_[0] * other.data_[0] + data_[1] * other.data_[1] + data_[2] * other.data_[2];
    }
    Vector operator*(double multiplier) const {
        return Vector(data_[0] * multiplier, data_[1] * multiplier, data_[2] * multiplier);
    }
    Vector FindColor(const Vector& other) const {
        return Vector{data_[0] * other.data_[0], data_[1] * other.data_[1],
                      data_[2] * other.data_[2]};
    }
    Vector operator+(const Vector& other) const {
        return Vector(data_[0] + other.data_[0], data_[1] + other.data_[1],
                      data_[2] + other.data_[2]);
    }
    Vector operator+(const double add) const {
        return Vector(data_[0] + add, data_[1] + add, data_[2] + add);
    }
    Vector operator-(const Vector& other) const {
        return *this + (other * (-1));
    }
    bool operator!=(const Vector& other) const {
        return data_[0] != other.data_[0] || data_[1] != other.data_[1] ||
               data_[2] != other.data_[2];
    }

    void Normalize() {
        double length = sqrt(*this * *this);

        if (length == 0) {
            throw std::runtime_error("Cannot normalize a zero-length vector.");
        }

        data_[0] /= length;
        data_[1] /= length;
        data_[2] /= length;
    }

private:
    std::array<double, 3> data_;
};

double DotProduct(const Vector& a, const Vector& b) {
    return a * b;
}
Vector CrossProduct(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
double Length(const Vector& v) {
    return sqrt(v * v);
}
