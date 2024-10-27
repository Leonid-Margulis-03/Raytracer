#include "test_asan.cpp"

TEST_CASE("Classic box") {
    CameraOptions camera_opts{.screen_width = 500,
                              .screen_height = 500,
                              .look_from = {-.5, 1.5, .98},
                              .look_to = {0., 1., 0.}};
    RenderOptions render_opts{4, RenderMode::kDepth};
    auto obj_filename = "classic_box/CornellBox.obj";

    std::filesystem::path output_path =
        "/Users/leo/Documents/YSDA/C++/margulisleo/raytracer-debug/output.png";
    CheckImage(obj_filename, "classic_box/depth1.png", camera_opts, render_opts, output_path);
    render_opts.mode = RenderMode::kNormal;
    CheckImage(obj_filename, "classic_box/normal1.png", camera_opts, render_opts, output_path);

    camera_opts.look_from = {-.9, 1.9, -1};
    camera_opts.look_to = {0., 0., 0.};
    render_opts.mode = RenderMode::kDepth;
    CheckImage(obj_filename, "classic_box/depth2.png", camera_opts, render_opts, output_path);
    render_opts.mode = RenderMode::kNormal;
    CheckImage(obj_filename, "classic_box/normal2.png", camera_opts, render_opts);
}

TEST_CASE("Deer") {
    CameraOptions camera_opts{.screen_width = 500,
                              .screen_height = 500,
                              .look_from = {100., 200., 150.},
                              .look_to = {0., 100., 0.}};
    RenderOptions render_opts{1, RenderMode::kDepth};
    std::filesystem::path output_path =
        "/Users/leo/Documents/YSDA/C++/margulisleo/raytracer-debug/output.png";
    CheckImage("deer/CERF_Free.obj", "deer/depth.png", camera_opts, render_opts, output_path);
    render_opts.mode = RenderMode::kNormal;
    CheckImage("deer/CERF_Free.obj", "deer/normal.png", camera_opts, render_opts, output_path);
}
