
#include "../include/gl.hpp"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

TEST_CASE( "rgb-hsl conversion", "UtilFunctions" ) {

    using Color = Eigen::Matrix<float, 3, 1>;

    // ---------------- DEG ----------------
    { //black
        const Color c = { 0., 0., 0. };
        const Color c_rgb = GL::HSL_to_RGB_deg(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 0);
    }
    { //white
        const Color c = { 0., 0., 1. };
        const Color c_rgb = GL::HSL_to_RGB_deg(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 255);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 255);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 255);
    }
    { //red
        const Color c = { 0., 1., 0.5 };
        const Color c_rgb = GL::HSL_to_RGB_deg(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 255);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 0);
    }
    { //green
        const Color c = { 120, 1., 0.5 };
        const Color c_rgb = GL::HSL_to_RGB_deg(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 255);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 0);
    }
    { //blue
        const Color c = { 240., 1., 0.5 };
        const Color c_rgb = GL::HSL_to_RGB_deg(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 255);
    }

    // ---------------- RAD ----------------
    { //black
        const Color c = { 0., 0., 0. };
        const Color c_rgb = GL::HSL_to_RGB_rad(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 0);
    }
    { //white
        const Color c = { 0., 0., 1. };
        const Color c_rgb = GL::HSL_to_RGB_rad(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 255);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 255);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 255);
    }
    { //red
        const Color c = { 0., 1., 0.5 };
        const Color c_rgb = GL::HSL_to_RGB_rad(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 255);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 0);
    }
    { //green
        const Color c = { 120 * (180. / M_PI), 1., 0.5 };
        const Color c_rgb = GL::HSL_to_RGB_rad(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 255);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 0);
    }
    { //blue
        const Color c = { 240. * (180. / M_PI), 1., 0.5 };
        const Color c_rgb = GL::HSL_to_RGB_rad(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 255);
    }
}

TEST_CASE( "complex to hsl conversion", "UtilFunctions" ) {

    using Color = Eigen::Matrix<float, 3, 1>;
    
    std::complex<float> cn (0.f, 1.f);
    Color col = GL::c_to_HSL(cn);
    const Color c_rgb = GL::HSL_to_RGB_deg(col);

    //REQUIRE(int(std::round(c_rgb(0) * 255)) == 0);

}
