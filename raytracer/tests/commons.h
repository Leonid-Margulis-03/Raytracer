#pragma once

#include "image.h"

#include <cmath>
#include <string>
#include <ranges>

#include <catch2/catch_test_macros.hpp>

double PixelDistance(const RGB& lhs, const RGB& rhs) {
    auto dr = lhs.r - rhs.r;
    auto dg = lhs.g - rhs.g;
    auto db = lhs.b - rhs.b;
    return std::sqrt(dr * dr + dg * dg + db * db);
}

void Compare(const Image& actual, const Image& expected) {
    constexpr auto kEps = 2.;
    auto matches = 0;

    REQUIRE(actual.Width() == expected.Width());
    REQUIRE(actual.Height() == expected.Height());
    for (auto y : std::views::iota(0, actual.Height())) {
        for (auto x : std::views::iota(0, actual.Width())) {
            auto actual_data = actual.GetPixel(y, x);
            auto expected_data = expected.GetPixel(y, x);
            auto diff = PixelDistance(actual_data, expected_data);
            matches += (diff < kEps);
        }
    }

    auto similarity = static_cast<double>(matches) / (actual.Width() * actual.Height());
    CHECK(similarity >= .99);
}
