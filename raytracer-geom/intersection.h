#pragma once

#include "vector.h"

class Intersection {
public:
    Intersection(const Vector& position, const Vector& normal, double distance)
        : position_(position), normal_(normal), distance_(distance) {
        this->normal_.Normalize();
    }
    const Vector& GetPosition() const {
        return position_;
    }
    const Vector& GetNormal() const {
        return normal_;
    }
    double GetDistance() const {
        return distance_;
    }

private:
    Vector position_;
    Vector normal_;
    double distance_;
};
