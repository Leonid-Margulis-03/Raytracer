#include "options/camera_options.h"
#include "options/render_options.h"
#include "tests/commons.h"
#include "raytracer.h"
#include "util.h"

#include <cmath>
#include <string_view>
#include <optional>

#include <catch2/catch_test_macros.hpp>

void CheckImage(std::string_view obj_filename, std::string_view result_filename,
                const CameraOptions& camera_options, const RenderOptions& render_options,
                const std::optional<std::filesystem::path>& output_path = std::nullopt) {
    static const auto kTestsDir = GetRelativeDir(__FILE__, "tests");
    auto image = Render(kTestsDir / obj_filename, camera_options, render_options);
    if (output_path) {
        image.Write(*output_path);
    }
    Compare(image, Image{kTestsDir / result_filename});
}

TEST_CASE("Shading parts") {
    CameraOptions camera_opts{640, 480};
    RenderOptions render_opts{1, RenderMode::kDepth};
    std::filesystem::path output_path =
        "/Users/leo/Documents/YSDA/C++/margulisleo/raytracer-debug/output.png";
    CheckImage("shading_parts/scene.obj", "shading_parts/depth.png", camera_opts, render_opts,
               output_path);
    render_opts.mode = RenderMode::kNormal;
    CheckImage("shading_parts/scene.obj", "shading_parts/normal.png", camera_opts, render_opts,
               output_path);
}

TEST_CASE("Triangle") {
    CameraOptions camera_opts{.screen_width = 640,
                              .screen_height = 480,
                              .look_from = {0., 2., 0.},
                              .look_to = {0., 0., 0.}};
    RenderOptions render_opts{1, RenderMode::kDepth};
    std::filesystem::path output_path =
        "/Users/leo/Documents/YSDA/C++/margulisleo/raytracer-debug/output.png";
    CheckImage("triangle/scene.obj", "triangle/depth.png", camera_opts, render_opts, output_path);
    render_opts.mode = RenderMode::kNormal;
    CheckImage("triangle/scene.obj", "triangle/normal.png", camera_opts, render_opts, output_path);
}

TEST_CASE("Triangle2") {
    CameraOptions camera_opts{.screen_width = 640,
                              .screen_height = 480,
                              .look_from = {0., -2., 0.},
                              .look_to = {0., 0., 0.}};
    RenderOptions render_opts{1, RenderMode::kNormal};
    CheckImage("triangle/scene.obj", "triangle/depth2.png", camera_opts, render_opts);
    render_opts.mode = RenderMode::kNormal;
    CheckImage("triangle/scene.obj", "triangle/normal2.png", camera_opts, render_opts);
}

TEST_CASE("Box with spheres") {
    CameraOptions camera_opts{.screen_width = 640,
                              .screen_height = 480,
                              .fov = std::numbers::pi / 3,
                              .look_from = {0., .7, 1.75},
                              .look_to = {0., .7, 0.}};
    RenderOptions render_opts{4, RenderMode::kDepth};
    std::filesystem::path output_path =
        "/Users/leo/Documents/YSDA/C++/margulisleo/raytracer-debug/output.png";
    CheckImage("box/cube.obj", "box/depth.png", camera_opts, render_opts, output_path);
    render_opts.mode = RenderMode::kNormal;
    CheckImage("box/cube.obj", "box/normal.png", camera_opts, render_opts, output_path);
}
