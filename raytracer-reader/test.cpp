#include "scene.h"
#include "util.h"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

auto WithinAbs(double x) {
    return Catch::Matchers::WithinAbs(x, 1e-12);
}

void Check(const Vector& vector, double x, double y, double z) {
    CHECK_THAT(vector[0], WithinAbs(x));
    CHECK_THAT(vector[1], WithinAbs(y));
    CHECK_THAT(vector[2], WithinAbs(z));
}

void Check(const Vector& vector, double x) {
    Check(vector, x, x, x);
}

TEST_CASE("Scene") {
    const auto test_dir = GetRelativeDir(__FILE__, "tests");
    const auto scene = ReadScene(test_dir / "cube.obj");

    const auto& materials_map = scene.GetMaterials();
    REQUIRE(materials_map.size() == 9);

    // objects
    const auto& objects = scene.GetObjects();
    REQUIRE(objects.size() == 10);

    {
        const auto& obj = objects[0];
        Check(obj.polygon[0], 1., 0., -1.04);
        Check(obj.polygon[1], -.99, 0., -1.04);
        Check(obj.polygon[2], -1.01, 0., .99);
        CHECK(obj.material->name == "floor");

        Check(*obj.GetNormal(0), .0, 1., 0.);
        Check(*obj.GetNormal(1), .0, 1., 0.);
        Check(*obj.GetNormal(2), .0, 1., 0.);
    }

    {
        const auto& obj = objects[4];
        Check(obj.polygon[0], 1., 1.59, -1.04);
        Check(obj.polygon[1], -1.02, 1.59, -1.04);
        Check(obj.polygon[2], -.99, 0., -1.04);
        CHECK(obj.material->name == "backWall");

        Check(*obj.GetNormal(0), 0., 0., 1.);
        Check(*obj.GetNormal(1), 0., 0., 1.);
        Check(*obj.GetNormal(2), 0., 0., 1.);
    }

    {
        const auto& obj = objects[5];
        Check(obj.polygon[0], 1., 1.59, -1.04);
        Check(obj.polygon[1], -.99, 0., -1.04);
        Check(obj.polygon[2], 1., 0., -1.04);
        CHECK(obj.material->name == "backWall");

        Check(*obj.GetNormal(0), 0., 0., 1.);
        Check(*obj.GetNormal(1), 0., 0., 1.);
        Check(*obj.GetNormal(2), 0., 0., 1.);
    }

    {
        const auto& obj = objects[8];
        Check(obj.polygon[0], -1.02, 1.59, -1.04);
        Check(obj.polygon[1], -1.02, 1.59, .99);
        Check(obj.polygon[2], -1.01, 0., .99);
        CHECK(obj.material->name == "leftWall");

        Check(*obj.GetNormal(0), .9999, .0135, .0057);
        Check(*obj.GetNormal(1), 1., .0063, 0.);
        Check(*obj.GetNormal(2), .9999, .0116, .0042);
    }

    Check(*objects[1].GetNormal(1), 0., 1., 0.);
    Check(*objects[6].GetNormal(2), -1., 0., 0.);

    for (const auto& object : objects) {
        CHECK(materials_map.contains(object.material->name));
    }

    // spheres
    const auto& spheres = scene.GetSphereObjects();
    REQUIRE(spheres.size() == 2);
    {
        const auto& obj = spheres[0];
        Check(obj.sphere.GetCenter(), -.4, .3, -.4);
        CHECK_THAT(obj.sphere.GetRadius(), WithinAbs(.3));
        CHECK(obj.material->name == "leftSphere");
    }
    {
        const auto& obj = spheres[1];
        Check(obj.sphere.GetCenter(), .3, .3, 0.);
        CHECK_THAT(obj.sphere.GetRadius(), WithinAbs(.3));
        CHECK(obj.material->name == "rightSphere");
    }

    // lights
    const auto& lights = scene.GetLights();
    REQUIRE(lights.size() == 2);
    Check(lights[0].position, 0., 1.5899, 0.);
    Check(lights[0].intensity, 1.);
    Check(lights[1].position, 0., .7, 1.98);
    Check(lights[1].intensity, .5);

    // materials
    const auto& right_sphere = materials_map.at("rightSphere");
    CHECK_THAT(right_sphere.specular_exponent, WithinAbs(1024.));
    CHECK_THAT(right_sphere.refraction_index, WithinAbs(1.8));
    Check(right_sphere.ambient_color, 0.);
    Check(right_sphere.diffuse_color, 0.);
    Check(right_sphere.specular_color, 0.);
    Check(right_sphere.albedo, 0., .3, .7);
    Check(right_sphere.intensity, 0.);

    const auto& light = materials_map.at("light");
    CHECK_THAT(light.specular_exponent, WithinAbs(10.));
    CHECK_THAT(light.refraction_index, WithinAbs(1.5));
    Check(light.ambient_color, .78);
    Check(light.diffuse_color, .78);
    Check(light.specular_color, 0.);
    Check(light.intensity, 1.);

    const auto& wall_behind = materials_map.at("wallBehind");
    CHECK_THAT(wall_behind.refraction_index, WithinAbs(1.));
    Check(wall_behind.diffuse_color, .2, .7, .8);
    Check(wall_behind.albedo, 1., 0., 0.);

    const auto& back_wall = materials_map.at("backWall");
    Check(back_wall.specular_color, 1.5);
    Check(back_wall.albedo, .5, 0., 0.);
    Check(back_wall.diffuse_color, .725, .91, .88);
}
