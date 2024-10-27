#include "geometry.h"
#include "vector.h"
#include "sphere.h"
#include "intersection.h"
#include "triangle.h"
#include "ray.h"
#include "util.h"

#include <cmath>
#include <optional>
#include <numbers>
#include <fstream>
#include <algorithm>
#include <array>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

namespace {

const auto kTestsDir = GetRelativeDir(__FILE__, "tests");

constexpr auto kX = 12.3;
constexpr auto kY = 45.6;
constexpr auto kZ = 78.9;

auto WithinAbs(double target, double eps = 1e-9) {
    return Catch::Matchers::WithinAbs(target, eps);
}

void CheckLength(Vector v, double expected) {
    CHECK_THAT(Length(v), WithinAbs(expected));
    v.Normalize();
    CHECK_THAT(Length(v), WithinAbs(1.));
}

Vector ReadVector(std::istream* is) {
    double x, y, z;
    *is >> x >> y >> z;
    return {x, y, z};
}

Sphere ReadSphere(std::istream* is) {
    auto center = ReadVector(is);
    double radius;
    *is >> radius;
    return {center, radius};
}

Triangle ReadTriangle(std::istream* is) {
    return {ReadVector(is), ReadVector(is), ReadVector(is)};
}

Ray ReadRay(std::istream* is) {
    return {ReadVector(is), ReadVector(is)};
}

void CheckEquals(const Vector& v, const std::array<double, 3>& expected) {
    CHECK(v[0] == expected[0]);
    CHECK(v[1] == expected[1]);
    CHECK(v[2] == expected[2]);
}

void CheckWithinAbs(const Vector& v, const Vector& expected) {
    CHECK_THAT(v[0], WithinAbs(expected[0]));
    CHECK_THAT(v[1], WithinAbs(expected[1]));
    CHECK_THAT(v[2], WithinAbs(expected[2]));
}

void CheckCoords(const Triangle& tr, const Vector& point, const std::array<double, 3>& coeffs) {
    std::array indexes = {0, 1, 2};
    do {
        auto [i0, i1, i2] = indexes;
        Triangle t{tr[i0], tr[i1], tr[i2]};
        Vector p{point[i0], point[i1], point[i2]};
        Vector expected = {coeffs[i0], coeffs[i1], coeffs[i2]};
        CheckWithinAbs(GetBarycentricCoords(t, point), expected);
    } while (std::ranges::next_permutation(indexes).found);
}

void CheckIntersection(const std::optional<Intersection>& intersection, std::istream* is) {
    double result;
    *is >> result;

    if ((result >= 0) != intersection.has_value()) {
        FAIL_CHECK("Wrong GetIntersection");
        return;
    }

    if (intersection) {
        CHECK_THAT(intersection->GetDistance(), WithinAbs(result));
        CheckWithinAbs(intersection->GetPosition(), ReadVector(is));
        CheckWithinAbs(intersection->GetNormal(), ReadVector(is));
    }
}

}  // namespace

TEST_CASE("Initialize vector") {
    CheckEquals({kX, kY, kZ}, {kX, kY, kZ});
    Vector vec;
    vec[0] = kX;
    vec[1] = kY;
    vec[2] = kZ;
    CheckEquals(vec, {kX, kY, kZ});
}

TEST_CASE("Normalize, Length") {
    {
        Vector vec{kX, 0, 0};
        CHECK_THAT(Length(vec), WithinAbs(kX));
        vec.Normalize();
        CheckWithinAbs(vec, {1, 0, 0});
    }
    {
        Vector vec{0, kY, 0};
        CHECK_THAT(Length(vec), WithinAbs(kY));
        vec.Normalize();
        CheckWithinAbs(vec, {0, 1, 0});
    }
    {
        Vector vec{0, 0, kZ};
        CHECK_THAT(Length(vec), WithinAbs(kZ));
        vec.Normalize();
        CheckWithinAbs(vec, {0, 0, 1});
    }

    CheckLength({3, 4, 0}, 5);
    CheckLength({56, 0, 33}, 65);
    CheckLength({0, 36, 77}, 85);
    CheckLength({10, 10, 23}, 27);
    CheckLength({3, 16, 24}, 29);
}

TEST_CASE("Dot product") {
    CHECK_THAT(DotProduct({kX, 1, 0}, {1, -kX, 0}), WithinAbs(0.));
    CHECK_THAT(DotProduct({0, kY, -1}, {0, 1, kY}), WithinAbs(0.));
    CHECK_THAT(DotProduct({3, kX, -kZ}, {kZ, kZ, 3}), WithinAbs(kX * kZ));

    Vector vec{kX, kY, kZ};
    CHECK_THAT(DotProduct(vec, {1, 0, 0}), WithinAbs(kX));
    CHECK_THAT(DotProduct(vec, {0, 2, 0}), WithinAbs(2 * kY));
    CHECK_THAT(DotProduct(vec, {0, 0, -3}), WithinAbs(-3 * kZ));
}

TEST_CASE("Cross product") {
    CheckWithinAbs(CrossProduct({kX, 0, 0}, {0, kY, 0}), {0, 0, kX * kY});
    CheckWithinAbs(CrossProduct({0, kY, 0}, {0, 0, kZ}), {kY * kZ, 0, 0});
    CheckWithinAbs(CrossProduct({0, 0, kZ}, {kX, 0, 0}), {0, kZ * kX, 0});
}

TEST_CASE("Triangle") {
    {
        Triangle triangle{{kX, 0, 0}, {0, kY, 0}, {0, 0, 0}};
        CheckEquals(triangle[0], {kX, 0, 0});
        CheckEquals(triangle[1], {0, kY, 0});
        CheckEquals(triangle[2], {0, 0, 0});
        CHECK_THAT(triangle.Area(), WithinAbs(.5 * kX * kY));
    }
    {
        Triangle triangle{{1, 0, 0}, {1, kY, 0}, {1, 0, kZ}};
        CheckEquals(triangle[0], {1, 0, 0});
        CheckEquals(triangle[1], {1, kY, 0});
        CheckEquals(triangle[2], {1, 0, kZ});
        CHECK_THAT(triangle.Area(), WithinAbs(.5 * kY * kZ));
    }
    {
        Triangle triangle{{0, 2, 1 + kZ}, {0, 2, 1}, {kX, 2, 1}};
        CheckEquals(triangle[0], {0, 2, 1 + kZ});
        CheckEquals(triangle[1], {0, 2, 1});
        CheckEquals(triangle[2], {kX, 2, 1});
        CHECK_THAT(triangle.Area(), WithinAbs(.5 * kZ * kX));
    }
}

TEST_CASE("Intersection") {
    Sphere sphere{{0, 0, 0}, 2};
    auto intersection = GetIntersection({{5, 0, 2.2}, {-1, 0, 0}}, sphere);
    CHECK_FALSE(intersection);

    intersection = GetIntersection({{5, 0, 0}, {-1, 0, 0}}, sphere);
    REQUIRE(intersection);
    CheckWithinAbs(intersection->GetPosition(), {2, 0, 0});
    CheckWithinAbs(intersection->GetNormal(), {1, 0, 0});
    CHECK_THAT(intersection->GetDistance(), WithinAbs(3.));

    intersection = GetIntersection({{5, 0, 2}, {-1, 0, 0}}, sphere);
    REQUIRE(intersection);
    CheckWithinAbs(intersection->GetPosition(), {0, 0, 2});
    CheckWithinAbs(intersection->GetNormal(), {0, 0, 1});
    CHECK_THAT(intersection->GetDistance(), WithinAbs(5.));

    intersection = GetIntersection({{0, 0, 0}, {-1, 0, 0}}, sphere);
    REQUIRE(intersection);
    CheckWithinAbs(intersection->GetPosition(), {-2, 0, 0});
    CheckWithinAbs(intersection->GetNormal(), {1, 0, 0});
    CHECK_THAT(intersection->GetDistance(), WithinAbs(2.));

    Triangle triangle{{0, 0, 0}, {4, 0, 0}, {0, 4, 0}};
    intersection = GetIntersection({{2, 2, 1}, {0, 0, -1}}, triangle);
    CheckWithinAbs(intersection->GetPosition(), {2, 2, 0});
    CheckWithinAbs(intersection->GetNormal(), {0, 0, 1});
    CHECK_THAT(intersection->GetDistance(), WithinAbs(1.));

    CHECK_FALSE(GetIntersection({{3, 3, 1}, {-1, -1, 0}}, triangle));
}

TEST_CASE("Sphere intersection") {
    std::ifstream is{kTestsDir / "sphere.txt"};
    int n;
    is >> n;

    while (n--) {
        auto ray = ReadRay(&is);
        auto sphere = ReadSphere(&is);
        CheckIntersection(GetIntersection(ray, sphere), &is);
    }
}

TEST_CASE("Triangle intersection") {
    std::ifstream is{kTestsDir / "triangle.txt"};
    int n;
    is >> n;

    while (n--) {
        auto ray = ReadRay(&is);
        auto triangle = ReadTriangle(&is);
        CheckIntersection(GetIntersection(ray, triangle), &is);
    }
}

TEST_CASE("Refract, Reflect") {
    Vector normal{0, 1, 0};
    auto d = std::numbers::sqrt2 / 2;
    Vector ray{d, -d, 0};
    CheckWithinAbs(Reflect(ray, normal), {d, d, 0});

    auto refract_opt = Refract(ray, normal, .9);
    REQUIRE(refract_opt);
    auto& refract = *refract_opt;

    CHECK_THAT(refract[0], WithinAbs(.636396, 1e-6));
    CHECK_THAT(refract[1], WithinAbs(-.771362, 1e-6));
    CHECK_THAT(refract[2], WithinAbs(0.));
}

TEST_CASE("Reflect") {
    std::ifstream is{kTestsDir / "reflect.txt"};
    int n;
    is >> n;

    while (n--) {
        auto v = ReadVector(&is);
        auto normal = ReadVector(&is);
        auto expected = ReadVector(&is);
        CheckWithinAbs(Reflect(v, normal), expected);
    }
}

TEST_CASE("Refract") {
    std::ifstream is{kTestsDir / "refract.txt"};
    int n;
    bool has_refract;
    is >> n;

    while (n--) {
        auto v = ReadVector(&is);
        auto normal = ReadVector(&is);
        auto result = Refract(v, normal, 1.1);
        is >> has_refract;
        if (result.has_value() != has_refract) {
            FAIL_CHECK("Wrong Refract");
        } else if (has_refract) {
            auto expected = ReadVector(&is);
            CheckWithinAbs(*result, expected);
        }
    }
}

TEST_CASE("Barycentric coords") {
    {
        Triangle t{{0, 0, 0}, {2, 0, 0}, {0, 2, 0}};
        CheckCoords(t, {1, 1, 0}, {0, .5, .5});
        CheckCoords(t, {2, 0, 0}, {0, 1, 0});
        CheckCoords(t, {.2, .2, 0}, {.8, .1, .1});
    }
    {
        Triangle t{{3, 0, 0}, {0, 3, 0}, {0, 0, 3}};
        CheckCoords(t, {0, 1.5, 1.5}, {0, .5, .5});
        CheckCoords(t, {0, 0, 3}, {0, 0, 1});
        CheckCoords(t, {1, 1, 1}, {1. / 3, 1. / 3, 1. / 3});
        CheckCoords(t, {2, 1, 0}, {2. / 3, 1. / 3, 0});
    }
    {
        Triangle t{{18, 4, 0}, {0, 15, 6}, {6, 0, 24}};
        CheckCoords(t, {8, 19. / 3, 10}, {1. / 3, 1. / 3, 1. / 3});
        CheckCoords(t, {10, 7, 6}, {3. / 6, 2. / 6, 1. / 6});
    }
}
