#pragma once

#include "../raytracer-geom/triangle.h"
#include "material.h"
#include "../raytracer-geom/sphere.h"
#include "../raytracer-geom/vector.h"

struct Object {
    Object(Material* material, Triangle polygon, Vector* normal1, Vector* normal2, Vector* normal3)
        : material(material),
          polygon(polygon),
          normal1_(normal1),
          normal2_(normal2),
          normal3_(normal3) {
    }
    const Material* material = nullptr;
    Triangle polygon;

    const Vector* GetNormal(size_t index) const {
        switch (index) {
            case 0:
                return normal1_;
            case 1:
                return normal2_;
            case 2:
                return normal3_;
            default:
                return nullptr;
        }
    }

private:
    Vector* normal1_ = nullptr;
    Vector* normal2_ = nullptr;
    Vector* normal3_ = nullptr;
};

struct SphereObject {
    const Material* material = nullptr;
    Sphere sphere;
};
