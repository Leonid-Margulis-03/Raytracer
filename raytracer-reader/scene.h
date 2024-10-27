#pragma once

#include "material.h"
#include "../raytracer-geom/vector.h"
#include "object.h"
#include "light.h"

#include <vector>
#include <unordered_map>
#include <string>
#include <filesystem>
#include <fstream>
#include <utility>
#include <sstream>

class Scene {
public:
    Scene(std::vector<Object> objects, std::vector<SphereObject> sphere_objects,
          std::vector<Light> lights, std::unordered_map<std::string, Material> materials,
          std::vector<Vector> vectors, std::vector<Vector> normals)
        : objects_(objects),
          sphere_objects_(sphere_objects),
          lights_(lights),
          materials_(materials),
          vectors_(vectors),
          normals_(normals) {
    }
    const std::vector<Object>& GetObjects() const {
        return this->objects_;
    }
    const std::vector<SphereObject>& GetSphereObjects() const {
        return this->sphere_objects_;
    }
    const std::vector<Light>& GetLights() const {
        return this->lights_;
    }
    const std::unordered_map<std::string, Material>& GetMaterials() const {
        return this->materials_;
    }

public:
    std::vector<Object> objects_;
    std::vector<SphereObject> sphere_objects_;
    std::vector<Light> lights_;
    std::unordered_map<std::string, Material> materials_;
    std::vector<Vector> vectors_;
    std::vector<Vector> normals_;
};

void ReadMaterials(const std::filesystem::path& path,
                   std::unordered_map<std::string, Material>& materials) {
    std::ifstream material_file(path);
    std::string line;
    std::string name;

    while (std::getline(material_file, line)) {
        // reading the command
        std::stringstream s(line);
        std::string command;
        s >> command;

        if (command == "newmtl") {
            s >> name;
            materials[name] = Material();
            materials[name].name = name;
            materials[name].albedo = Vector(1, 0, 0);
            materials[name].refraction_index = 1;
        }

        if (command == "Ka") {
            double r, g, b;
            s >> r >> g >> b;
            materials[name].ambient_color = Vector(r, g, b);
        }

        if (command == "Kd") {
            double r, g, b;
            s >> r >> g >> b;
            materials[name].diffuse_color = Vector(r, g, b);
        }

        if (command == "Ks") {
            double r, g, b;
            s >> r >> g >> b;
            materials[name].specular_color = Vector(r, g, b);
        }

        if (command == "Ke") {
            double r, g, b;
            s >> r >> g >> b;
            materials[name].intensity = Vector(r, g, b);
        }

        if (command == "Ns") {
            double specular_exp;
            s >> specular_exp;
            materials[name].specular_exponent = specular_exp;
        }

        if (command == "Ni") {
            double refract_index;
            s >> refract_index;
            materials[name].refraction_index = refract_index;
        }

        if (command == "al") {
            double a, b, c;
            s >> a >> b >> c;
            materials[name].albedo = Vector(a, b, c);
        }
    }
}
Scene ReadScene(const std::filesystem::path& path) {
    std::ifstream obj_file(path);
    std::string line;
    std::filesystem::path obj_directory = path.parent_path();

    Scene scene{std::vector<Object>(), std::vector<SphereObject>(),
                std::vector<Light>(),  std::unordered_map<std::string, Material>(),
                std::vector<Vector>(), std::vector<Vector>()};
    scene.vectors_.reserve(1000);
    scene.normals_.reserve(1000);
    scene.objects_.reserve(1000);
    scene.sphere_objects_.reserve(1000);
    scene.lights_.reserve(1000);

    Material* material = nullptr;
    while (std::getline(obj_file, line)) {
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::stringstream s(line);
        std::string command;
        s >> command;
        // Include new materials
        if (command == "mtllib") {
            std::string path_mtl;
            s >> path_mtl;
            std::filesystem::path path_material = obj_directory / path_mtl;
            ReadMaterials(path_material, scene.materials_);
        }
        // Vector
        if (command == "v") {
            double x, y, z;
            s >> x >> y >> z;
            scene.vectors_.push_back(Vector{x, y, z});
        }
        // Vector normal
        if (command == "vn") {
            double x, y, z;
            s >> x >> y >> z;
            scene.normals_.push_back(Vector{x, y, z});
            //  scene.normals_.back().Normalize();
        }

        if (command == "f") {
            std::string number;
            int count = 0;

            int v1 = 0;
            int vn1 = 0;
            int v2 = 0;
            int vn2 = 0;
            int v3 = 0;
            int vn3 = 0;
            int count_slash = 0;
            for (int i = 2; i < static_cast<int>(line.length()); ++i) {
                if (line[i] == '/') {
                    count_slash += 1;
                }
            }
            if (count_slash == 0) {
                s >> v1 >> v2 >> v3;

                if (v1 < 0) {
                    v1 = static_cast<int>(scene.vectors_.size()) + v1 + 1;
                }
                if (v2 < 0) {
                    v2 = static_cast<int>(scene.vectors_.size()) + v2 + 1;
                }
                if (v3 < 0) {
                    v3 = static_cast<int>(scene.vectors_.size()) + v3 + 1;
                }

                scene.objects_.push_back(
                    Object{material,
                           Triangle{scene.vectors_[v1 - 1], scene.vectors_[v2 - 1],
                                    scene.vectors_[v3 - 1]},
                           nullptr, nullptr, nullptr});
                int tmp = 0;
                while (s >> tmp) {
                    v2 = v3;
                    v3 = tmp;
                    if (v3 < 0) {
                        v3 = static_cast<int>(scene.vectors_.size()) + v3 + 1;
                    }
                    scene.objects_.push_back(
                        Object{material,
                               Triangle{scene.vectors_[v1 - 1], scene.vectors_[v2 - 1],
                                        scene.vectors_[v3 - 1]},
                               nullptr, nullptr, nullptr});
                }
            } else {
                for (int i = 2; i < static_cast<int>(line.length()); ++i) {
                    if (line[i] == ' ') {
                        if (count == 2) {
                            if (vn1 == 0) {
                                vn1 = stoi(number);
                                if (vn1 < 0) {
                                    vn1 = static_cast<int>(scene.normals_.size()) + vn1 + 1;
                                }
                            } else {
                                vn2 = vn3;
                                vn3 = stoi(number);
                                if (vn3 < 0) {
                                    vn3 = static_cast<int>(scene.normals_.size()) + vn3 + 1;
                                }
                            }
                        }
                        number = "";
                        // if we have 3 vectors
                        if (v2 != 0) {

                            if (vn1 != 0) {
                                scene.objects_.push_back(
                                    Object{material,
                                           Triangle{scene.vectors_[v1 - 1], scene.vectors_[v2 - 1],
                                                    scene.vectors_[v3 - 1]},
                                           &scene.normals_[vn1 - 1], &scene.normals_[vn2 - 1],
                                           &scene.normals_[vn3 - 1]});
                            } else {
                                scene.objects_.push_back(
                                    Object{material,
                                           Triangle{scene.vectors_[v1 - 1], scene.vectors_[v2 - 1],
                                                    scene.vectors_[v3 - 1]},
                                           nullptr, nullptr, nullptr});
                            }
                        }
                        count = 0;
                        continue;
                    }
                    if (line[i] == '/') {
                        if (count == 0) {
                            if (v1 == 0) {
                                v1 = stoi(number);
                                if (v1 < 0) {
                                    v1 = static_cast<int>(scene.vectors_.size()) + v1 + 1;
                                }
                            } else {
                                v2 = v3;
                                v3 = stoi(number);
                                if (v3 < 0) {
                                    v3 = static_cast<int>(scene.vectors_.size()) + v3 + 1;
                                }
                            }
                        }
                        number = "";
                        count++;
                        continue;
                    }
                    number += line[i];
                }
                if (count == 2) {
                    if (vn1 == 0) {
                        vn1 = stoi(number);
                        if (vn1 < 0) {
                            vn1 = static_cast<int>(scene.normals_.size()) + vn1 + 1;
                        }
                    } else {
                        vn2 = vn3;
                        vn3 = stoi(number);
                        if (vn3 < 0) {
                            vn3 = static_cast<int>(scene.normals_.size()) + vn3 + 1;
                        }
                    }
                }
                number = "";
                if (count != 0) {
                    if (vn1 != 0) {
                        scene.objects_.push_back(
                            Object{material,
                                   Triangle{scene.vectors_[v1 - 1], scene.vectors_[v2 - 1],
                                            scene.vectors_[v3 - 1]},
                                   &scene.normals_[vn1 - 1], &scene.normals_[vn2 - 1],
                                   &scene.normals_[vn3 - 1]});
                    } else {
                        scene.objects_.push_back(
                            Object{material,
                                   Triangle{scene.vectors_[v1 - 1], scene.vectors_[v2 - 1],
                                            scene.vectors_[v3 - 1]},
                                   nullptr, nullptr, nullptr});
                    }
                }
            }
        }
        // Change of material
        if (command == "usemtl") {
            std::string material_name;
            s >> material_name;
            material = &scene.materials_[material_name];
        }
        // Sphere
        if (command == "S") {
            double x, y, z, r;
            s >> x >> y >> z >> r;

            scene.sphere_objects_.push_back(SphereObject{material, Sphere{Vector{x, y, z}, r}});
        }
        // Light
        if (command == "P") {
            double x, y, z, r, g, b;
            s >> x >> y >> z >> r >> g >> b;
            scene.lights_.push_back(Light{Vector{x, y, z}, Vector{r, g, b}});
        }
    }

    return scene;
}
