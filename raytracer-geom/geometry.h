#pragma once

#include "vector.h"
#include "sphere.h"
#include "intersection.h"
#include "triangle.h"
#include "ray.h"
#include <iostream>

#include <optional>

double Determinant(const Vector& v1, const Vector& v2, const Vector& v3) {
    return v1[0] * (v2[1] * v3[2] - v2[2] * v3[1]) + v2[0] * (v3[1] * v1[2] - v3[2] * v1[1]) +
           v3[0] * (v1[1] * v2[2] - v2[1] * v1[2]);
}
std::optional<Intersection> GetIntersection(const Ray& ray, const Sphere& sphere) {
    double r = sphere.GetRadius();
    Vector oc = ray.GetOrigin() - sphere.GetCenter();
    double a = ray.GetDirection() * ray.GetDirection();
    double b = 2 * (oc * ray.GetDirection());
    double c = oc * oc - r * r;

    double discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return std::nullopt;
    }
    double t1 = (-b - sqrt(discriminant)) / (2 * a);
    double t2 = (-b + sqrt(discriminant)) / (2 * a);

    double ans = std::max(t1, t2);

    if (t1 >= 0 && t1 < ans) {
        ans = t1;
    }
    if (t2 >= 0 && t2 < ans) {
        ans = t2;
    }

    if (ans < 0) {
        return std::nullopt;
    } else {
        Vector intersection = ray.GetOrigin() + ray.GetDirection() * ans;
        Vector normal = intersection - sphere.GetCenter();

        if (normal * ray.GetDirection() > 0) {
            normal = normal * (-1);
        }

        normal.Normalize();
        double distance = Length(intersection - ray.GetOrigin());
        return Intersection{intersection, normal, distance};
    }
}
std::optional<Intersection> GetIntersection(const Ray& ray, const Triangle& triangle) {
    // check for parallel
    Vector normal = CrossProduct(triangle[1] - triangle[0], triangle[2] - triangle[0]);
    normal.Normalize();

    if (normal * ray.GetDirection() > 0) {
        normal = normal * (-1);
    }

    if (ray.GetDirection() * normal == 0) {
        return std::nullopt;
    } else {
        Vector b = ray.GetOrigin() - triangle[0];
        double determinant = Determinant(ray.GetDirection() * (-1), triangle[1] - triangle[0],
                                         triangle[2] - triangle[0]);
        double t =
            Determinant(b, triangle[1] - triangle[0], triangle[2] - triangle[0]) / determinant;
        double u =
            Determinant(ray.GetDirection() * (-1), b, triangle[2] - triangle[0]) / determinant;
        double v =
            Determinant(ray.GetDirection() * (-1), triangle[1] - triangle[0], b) / determinant;
        if (t >= 0 && u >= 0 && v >= 0 && u + v <= 1) {
            double w = 1 - u - v;
            Vector intersection = triangle[0] * w + triangle[1] * u + triangle[2] * v;
            return Intersection{intersection, normal, Length(ray.GetOrigin() - intersection)};
        } else {
            return std::nullopt;
        }
    }
}

Vector Reflect(const Vector& ray, const Vector& normal) {
    double dot_product = -1 * (ray * normal);
    return ray + normal * dot_product * 2;
}
std::optional<Vector> Refract(const Vector& ray, const Vector& normal, double eta) {
    double cos_alpha = -1 * (ray * normal);
    double sin_alpha = sqrt(1 - cos_alpha * cos_alpha);
    if (sin_alpha * eta > 1) {
        return std::nullopt;
    }
    Vector result =
        ray * eta + normal * (eta * cos_alpha - sqrt(1 - eta * eta * sin_alpha * sin_alpha));
    return result;
}
Vector GetBarycentricCoords(const Triangle& triangle, const Vector& point) {
    // p = A*w + B*u + C*v
    double abc = triangle.Area();
    double pbc = Length(CrossProduct((triangle[1] - point), (triangle[2] - point))) * 1 / 2;
    double pca = Length(CrossProduct((triangle[2] - point), (triangle[0] - point))) * 1 / 2;

    double x = pbc / abc;
    double y = pca / abc;
    return Vector{x, y, 1 - x - y};
}
