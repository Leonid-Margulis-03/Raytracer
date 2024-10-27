#pragma once

#include "options/camera_options.h"
#include "options/render_options.h"
#include "image.h"
#include "scene.h"
#include "ray.h"
#include "geometry.h"
#include "material.h"
#include <cmath>
#include <filesystem>
#include <limits>
#include <iostream>

class Matrix {
public:
    Matrix(std::size_t rows, std::size_t columns)
        : matrix_(std::vector<std::vector<double>>(rows, std::vector<double>(columns, 0))),
          rows_(rows),
          columns_(columns) {
    }
    Matrix(std::size_t n)
        : matrix_(std::vector<std::vector<double>>(n, std::vector<double>(n, 0))),
          rows_(n),
          columns_(n) {
    }
    Matrix(const std::vector<std::vector<double>>& array)
        : matrix_(array), rows_(array.size()), columns_(array[0].size()) {
        for (int i = 0; i < static_cast<int>(array.size()); ++i) {
            for (int j = 0; j < static_cast<int>(array[0].size()); ++j) {
                (*this)(i, j) = array[i][j];
            }
        }
    }

    std::size_t Rows() const {
        return rows_;
    }

    std::size_t Columns() const {
        return columns_;
    }

    double& operator()(int row, int column) {
        return matrix_[row][column];
    }

    double operator()(int row, int column) const {
        return matrix_[row][column];
    }

    Matrix& operator+=(const Matrix& other) {
        for (std::size_t i = 0; i < other.rows_; ++i) {
            for (std::size_t j = 0; j < other.columns_; ++j) {
                this->matrix_[i][j] += other.matrix_[i][j];
            }
        }
        return *this;
    }

    Matrix& operator-=(const Matrix& other) {
        for (std::size_t i = 0; i < other.rows_; ++i) {
            for (std::size_t j = 0; j < other.columns_; ++j) {
                this->matrix_[i][j] -= other.matrix_[i][j];
            }
        }
        return *this;
    }

    Matrix operator+(const Matrix& other) const {
        Matrix tmp = *this;
        return tmp += other;
    }

    Matrix operator-(const Matrix& other) const {
        Matrix tmp = *this;
        return tmp -= other;
    }

    Matrix operator*(const Matrix& other) const {
        Matrix result{this->rows_, other.columns_};
        for (std::size_t i = 0; i < this->rows_; ++i) {
            for (std::size_t c = 0; c < other.columns_; ++c) {
                double sum = 0;
                for (std::size_t j = 0; j < this->columns_; ++j) {
                    sum += (*this)(i, j) * other(j, c);
                }
                result(i, c) = sum;
            }
        }
        return result;
    }

    Matrix& operator*=(const Matrix& other) {
        (*this) = (*this) * other;
        return *this;
    }

    Matrix operator*(const double m) const {
        Matrix result{this->rows_, this->columns_};
        for (std::size_t i = 0; i < this->rows_; ++i) {
            for (std::size_t c = 0; c < this->columns_; ++c) {
                result.matrix_[i][c] = matrix_[i][c] * m;
            }
        }
        return result;
    }
    Vector operator*(const Vector& v) const {
        Vector result = Vector{0, 0, 0};
        for (std::size_t i = 0; i < this->rows_; ++i) {
            double sum = 0;
            for (std::size_t c = 0; c < this->columns_; ++c) {
                sum += matrix_[i][c] * v[c];
            }
            result[i] = sum;
        }
        return result;
    }

private:
    std::vector<std::vector<double>> matrix_;
    std::size_t rows_;
    std::size_t columns_;
};

bool ApproximatelyEqual(double a, double b, double epsilon) {
    return fabs(a - b) < epsilon;
}
Matrix RotationMatrix(Vector& forward) {
    Vector z = forward;
    z.Normalize();
    Vector tmp = {0, 1, 0};
    if (z[0] == 0 && z[2] == 0 && z[1] == 1) {
        tmp = {0, 0, -1};
    }
    if (z[0] == 0 && z[2] == 0 && z[1] == -1) {
        tmp = {0, 0, 1};
    }
    Vector x = CrossProduct(tmp, z);
    x.Normalize();
    Vector y = CrossProduct(z, x);
    y.Normalize();

    if (DotProduct((CrossProduct(x, y)), z) < 0) {
        y = y * (-1);
    }
    assert(DotProduct((CrossProduct(x, y)), z) > 0);

    std::vector<std::vector<double>> matrix(3, std::vector<double>(3, 0));
    matrix[0][0] = x[0];
    matrix[1][0] = x[1];
    matrix[2][0] = x[2];

    matrix[0][1] = y[0];
    matrix[1][1] = y[1];
    matrix[2][1] = y[2];

    matrix[0][2] = z[0];
    matrix[1][2] = z[1];
    matrix[2][2] = z[2];
    Matrix rotation_matrix = Matrix(matrix);

    return rotation_matrix;
}

Vector FindRgbDebug(const Scene& scene, const Ray& ray_direction, int depth, int max_depth) {
    if (depth >= max_depth && max_depth != -1) {
        return Vector{0, 0, 0};
    }
    double closest_distance = std::numeric_limits<double>::max();
    Material material;
    Intersection intersection_point =
        Intersection(ray_direction.GetOrigin(), ray_direction.GetDirection(), -1);
    const Vector* normal_intersection_polygon = nullptr;
    // 2. intersections with polygons
    for (auto x : scene.objects_) {
        auto intersection = GetIntersection(ray_direction, x.polygon);

        if (intersection.has_value()) {
            Intersection tmp_intersection_point = intersection.value();
            if (tmp_intersection_point.GetDistance() < closest_distance) {
                normal_intersection_polygon = x.GetNormal(0);

                material = *x.material;
                intersection_point = intersection.value();
                closest_distance = tmp_intersection_point.GetDistance();
            }
        }
    }

    // 2. intersections with spheres
    for (auto x : scene.sphere_objects_) {
        auto intersection = GetIntersection(ray_direction, x.sphere);

        if (intersection.has_value()) {
            Intersection tmp_intersection_point = intersection.value();
            if (tmp_intersection_point.GetDistance() < closest_distance) {
                normal_intersection_polygon = nullptr;
                material = *x.material;
                intersection_point = intersection.value();
                closest_distance = intersection_point.GetDistance();
            }
        }
    }
    // then start recursion with depth
    /*
        std::string name;
        Vector ambient_color;
        Vector diffuse_color;
        Vector specular_color;
        Vector intensity;
        double specular_exponent;
        double refraction_index;
        Vector albedo;
     */
    // 3.          Ka                       Ke
    //    Vector Ibase = material.ambient_color + material.intensity;
    if (max_depth != -1) {
        double d = intersection_point.GetDistance();
        if (d == -1) {
            return Vector{1, 1, 1};
        } else {
            return Vector{d, d, d};
        }
    } else {
        // Normal mode
        double d = intersection_point.GetDistance();
        if (d == -1) {
            return Vector{0, 0, 0};
        } else {
            Vector normal = intersection_point.GetNormal();
            if (normal_intersection_polygon != nullptr) {
                normal = *normal_intersection_polygon;
            }
            normal = normal * 0.5;
            normal = normal + (0.5);
            return normal;
        }
    }

    //    Vector reflection = Reflect(ray_direction.GetDirection(), intersection_point.GetNormal());
    //    auto refraction = Refract(ray_direction.GetDirection(), intersection_point.GetNormal(),
    //    m.refraction_index);
}

Image RenderDepthNormal(const Scene& scene, const CameraOptions& camera_options, int depth) {
    int width_pixel = camera_options.screen_width;
    int height_pixel = camera_options.screen_height;

    Image image = Image(width_pixel, height_pixel);

    double height_screen = 2 * tan(camera_options.fov / 2);
    double width_screen = height_screen * width_pixel / height_pixel;

    double pixel_width = width_screen / width_pixel;
    double pixel_height = height_screen / height_pixel;
    double max_d = 0;
    std::vector<std::vector<Vector>> pixel_colors(width_pixel, std::vector<Vector>(height_pixel));
    Vector default_value = Vector(1, 1, 1);

    // 0. Find New Vector we look to
    Vector forward = (camera_options.look_from - camera_options.look_to);
    forward.Normalize();
    // 0.5 Find Rotation Matrix
    Matrix rotation_matrix = RotationMatrix(forward);
    for (int i = 0; i < width_pixel; ++i) {
        for (int j = 0; j < height_pixel; ++j) {
            double x = (i + 0.5) * pixel_width - (width_screen / 2);
            double y = -((j + 0.5) * pixel_height - (height_screen / 2));
            double z = -1;

            // 1. Camera CO
            Vector vector_direction = Vector(x, y, z);
            vector_direction.Normalize();
            // 2. Rotate our vector
            vector_direction = rotation_matrix * vector_direction;
            vector_direction.Normalize();
            // 3. Call RenderDepth with new rotated vector
            Vector rgb =
                FindRgbDebug(scene, Ray(camera_options.look_from, vector_direction), 0, depth);
            // 4. Save this vector to pixel_colors
            pixel_colors[i][j] = rgb;
            // 5. Find max distance to normalize FOR DEPTH
            if (pixel_colors[i][j] != default_value) {
                max_d = std::max(max_d, rgb[0]);
            }
        }
    }

    // 8. Normalization
    for (int i = 0; i < width_pixel; ++i) {
        for (int j = 0; j < height_pixel; ++j) {
            if (pixel_colors[i][j] != default_value && depth != -1) {
                pixel_colors[i][j] = pixel_colors[i][j] * (1 / max_d);
            }
            RGB rgb = {static_cast<int>(pixel_colors[i][j][0] * 255),
                       static_cast<int>(pixel_colors[i][j][1] * 255),
                       static_cast<int>(pixel_colors[i][j][2] * 255)};
            image.SetPixel(rgb, j, i);
        }
    }

    return image;
}

void ToneMappingGammaCorrection(std::vector<std::vector<Vector>>& pixel_colors) {
    double c = 0;
    // 0. find max C
    for (int i = 0; i < static_cast<int>(pixel_colors.size()); ++i) {
        for (int j = 0; j < static_cast<int>(pixel_colors[0].size()); ++j) {
            for (int k = 0; k < 3; ++k) {
                if (c < pixel_colors[i][j][k]) {
                    c = pixel_colors[i][j][k];
                }
            }
        }
    }
    // 1. Tone Mapping and Gamma Correction
    if (c != 0) {
        for (int i = 0; i < static_cast<int>(pixel_colors.size()); ++i) {
            for (int j = 0; j < static_cast<int>(pixel_colors[0].size()); ++j) {
                for (int k = 0; k < 3; ++k) {
                    double v_in = pixel_colors[i][j][k];
                    // Tone Mapping
                    double v_out = v_in * (1 + v_in / (c * c)) / (1 + v_in);
                    // Gamma Correction
                    double v_gamma = std::pow(v_out, 1.0 / 2.2);
                    pixel_colors[i][j][k] = v_gamma;
                }
            }
        }
    }
}

void FindIntersection(const Scene& scene, const Ray& ray_direction, Material& material,
                      const Vector*& normal_intersection_polygon, Intersection& intersection_point,
                      bool& flag) {
    double closest_distance = std::numeric_limits<double>::max();
    // 2. intersections with polygons
    for (auto x : scene.objects_) {
        auto intersection = GetIntersection(ray_direction, x.polygon);

        if (intersection.has_value()) {
            Intersection tmp_intersection_point = intersection.value();
            if (tmp_intersection_point.GetDistance() < closest_distance) {
                normal_intersection_polygon = x.GetNormal(0);

                material = *x.material;
                intersection_point = intersection.value();
                closest_distance = tmp_intersection_point.GetDistance();
            }
        }
    }

    // 2. intersections with spheres
    for (auto x : scene.sphere_objects_) {
        auto intersection = GetIntersection(ray_direction, x.sphere);

        if (intersection.has_value()) {
            Intersection tmp_intersection_point = intersection.value();
            if (tmp_intersection_point.GetDistance() < closest_distance) {
                normal_intersection_polygon = nullptr;
                flag = true;
                material = *x.material;
                intersection_point = intersection.value();
                closest_distance = intersection_point.GetDistance();
            }
        }
    }
}
Vector FindRgbFull(const Scene& scene, const Ray& ray_direction, int depth, bool flag_inside) {
    // If depth of recursion is achieved
    if (depth <= 0) {
        return Vector{0, 0, 0};
    }
    Material material;
    Intersection intersection_point =
        Intersection(ray_direction.GetOrigin(), ray_direction.GetDirection(), -1);
    const Vector* normal_intersection_polygon = nullptr;
    bool is_sphere = false;
    FindIntersection(scene, ray_direction, material, normal_intersection_polygon,
                     intersection_point, is_sphere);

    double distance = intersection_point.GetDistance();
    // if no intersection
    if (distance == -1) {
        return Vector{0, 0, 0};
    } else {
        Vector normal = intersection_point.GetNormal();
        if (normal_intersection_polygon != nullptr) {
            normal = *normal_intersection_polygon;
        }
        // then start recursion with depth
        /*
         *  Material
            std::string name;
            Vector ambient_color; Ka
            Vector diffuse_color; Kd
            Vector specular_color; Ks
            Vector intensity; Ke
            double specular_exponent; Ns
            double refraction_index; Ni
            Vector albedo;

            Light
            Vector position;
            Vector intensity;
         */
        Vector sum = {0, 0, 0};
        for (int i = 0; i < static_cast<int>(scene.lights_.size()); ++i) {
            Light l = scene.lights_[i];
            Vector l_m = l.position - intersection_point.GetPosition();
            l_m.Normalize();
            // start a ray from Light to current Point
            Ray ray_direction_from_light = Ray(l.position, l_m * (-1));
            // Find first intersection from this light towards the "point"
            Intersection intersection_point2 = Intersection(l.position, l_m * (-1), -1);

            Material material2;
            const Vector* normal_intersection_polygon2 = nullptr;
            FindIntersection(scene, ray_direction_from_light, material2,
                             normal_intersection_polygon2, intersection_point2, is_sphere);

            Vector normal2 = intersection_point2.GetNormal();
            if (normal_intersection_polygon2 != nullptr) {
                normal2 = *normal_intersection_polygon2;
            }
            normal2.Normalize();

            if (ApproximatelyEqual(intersection_point2.GetDistance(),
                                   Length(l.position - intersection_point.GetPosition()), 1e-9) &&
                DotProduct(intersection_point.GetNormal(), l_m) > 0) {
                double diffuse_intensity = DotProduct(l_m, normal2);
                if (diffuse_intensity > 0) {
                    sum = sum + material.diffuse_color.FindColor(l.intensity) *
                                    diffuse_intensity;  // light.intensity represents i_{m,d}
                }
                Vector r_m = Reflect(l_m * (-1), normal2);
                r_m.Normalize();
                double specular_intensity =
                    pow(std::max(0.0, DotProduct(r_m, ray_direction.GetDirection() * -1)),
                        material.specular_exponent);

                if (specular_intensity > 0) {
                    sum = sum + material2.specular_color.FindColor(l.intensity) *
                                    specular_intensity;  // light.intensity represents i_{m,s}
                }
            }
        }

        Vector al = material.albedo;
        Vector i_base = material.ambient_color + material.intensity + sum * al[0];
        Vector i_comp = Vector{0, 0, 0};

        // Reflected part
        Ray reflected = Ray(intersection_point.GetPosition() + normal * 1e-6,
                            Reflect(ray_direction.GetDirection(), normal));

        if (al[1] != 0 && !flag_inside) {
            i_comp = i_comp + FindRgbFull(scene, reflected, depth - 1, flag_inside) * al[1];
        }
        double t = al[2];
        if (flag_inside) {
            t = 1;
        }
        // Refracted part
        if (t != 0) {
            std::optional<Vector> vector_refracted = std::nullopt;
            if (!flag_inside) {
                vector_refracted =
                    Refract(ray_direction.GetDirection(), normal, 1 / material.refraction_index);
            }
            if (flag_inside) {
                vector_refracted =
                    Refract(ray_direction.GetDirection(), normal, material.refraction_index);
            }

            if (vector_refracted.has_value() && DotProduct(vector_refracted.value(), normal) < 0) {
                Ray refracted =
                    Ray(intersection_point.GetPosition() - normal * 1e-6, vector_refracted.value());
                if (is_sphere) {
                    i_comp = i_comp + FindRgbFull(scene, refracted, depth - 1, !flag_inside) * t;
                } else {
                    i_comp = i_comp + FindRgbFull(scene, refracted, depth - 1, flag_inside) * t;
                }
            }
        }
        return i_base + i_comp;
    }
}

Image RenderFull(const Scene& scene, const CameraOptions& camera_options, int depth) {
    int width_pixel = camera_options.screen_width;
    int height_pixel = camera_options.screen_height;

    Image image = Image(width_pixel, height_pixel);

    double height_screen = 2 * tan(camera_options.fov / 2);
    double width_screen = height_screen * width_pixel / height_pixel;

    double pixel_width = width_screen / width_pixel;
    double pixel_height = height_screen / height_pixel;

    std::vector<std::vector<Vector>> pixel_colors(width_pixel, std::vector<Vector>(height_pixel));

    // 0. Find New Vector we look to
    Vector forward = (camera_options.look_from - camera_options.look_to);
    forward.Normalize();
    // 0.5 Find Rotation Matrix
    Matrix rotation_matrix = RotationMatrix(forward);
    for (int i = 0; i < width_pixel; ++i) {
        for (int j = 0; j < height_pixel; ++j) {
            double x = (i + 0.5) * pixel_width - (width_screen / 2);
            double y = -((j + 0.5) * pixel_height - (height_screen / 2));
            double z = -1;

            // 1. Camera CO
            Vector vector_direction = Vector(x, y, z);
            vector_direction.Normalize();
            // 2. Rotate our vector
            vector_direction = rotation_matrix * vector_direction;
            vector_direction.Normalize();
            // 3. Call RenderDepth with new rotated vector
            Vector rgb =
                FindRgbFull(scene, Ray(camera_options.look_from, vector_direction), depth, false);
            // 4. Save this vector to pixel_colors
            pixel_colors[i][j] = rgb;
        }
    }
    // 6. Tone Mapping and GammaCorrection
    ToneMappingGammaCorrection(pixel_colors);

    // 7. Normalization to 255
    for (int i = 0; i < width_pixel; ++i) {
        for (int j = 0; j < height_pixel; ++j) {
            RGB rgb = {static_cast<int>(pixel_colors[i][j][0] * 255),
                       static_cast<int>(pixel_colors[i][j][1] * 255),
                       static_cast<int>(pixel_colors[i][j][2] * 255)};
            image.SetPixel(rgb, j, i);
        }
    }

    return image;
}

Image Render(const std::filesystem::path& path, const CameraOptions& camera_options,
             const RenderOptions& render_options) {
    Scene scene = ReadScene(path);  // build the scene
    if (render_options.mode == RenderMode::kDepth) {
        return RenderDepthNormal(scene, camera_options, render_options.depth);
    } else if (render_options.mode == RenderMode::kNormal) {
        return RenderDepthNormal(scene, camera_options, -1);
    } else {
        return RenderFull(scene, camera_options, render_options.depth);
    }
}
