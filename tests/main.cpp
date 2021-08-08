
#include "../include/hagedornrenderer.hpp"
#include "../include/math.hpp"
#include "../include/hdf5.h"

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

TEST_CASE( "rgb-hsl conversion", "UtilFunctions" ) {

    using Color = Eigen::Matrix<float, 3, 1>;

    // ---------------- DEG ----------------
    { //black
        const Color c = { 0., 0., 0. };
        const Color c_rgb = GL::Detail::HSL_to_RGB_deg(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 0);
    }
    { //white
        const Color c = { 0., 0., 1. };
        const Color c_rgb = GL::Detail::HSL_to_RGB_deg(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 255);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 255);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 255);
    }
    { //red
        const Color c = { 0., 1., 0.5 };
        const Color c_rgb = GL::Detail::HSL_to_RGB_deg(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 255);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 0);
    }
    { //green
        const Color c = { 120, 1., 0.5 };
        const Color c_rgb = GL::Detail::HSL_to_RGB_deg(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 255);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 0);
    }
    { //blue
        const Color c = { 240., 1., 0.5 };
        const Color c_rgb = GL::Detail::HSL_to_RGB_deg(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 255);
    }

    // ---------------- RAD ----------------
    { //black
        const Color c = { 0., 0., 0. };
        const Color c_rgb = GL::Detail::HSL_to_RGB_rad(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 0);
    }
    { //white
        const Color c = { 0., 0., 1. };
        const Color c_rgb = GL::Detail::HSL_to_RGB_rad(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 255);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 255);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 255);
    }
    { //red
        const Color c = { 0., 1., 0.5 };
        const Color c_rgb = GL::Detail::HSL_to_RGB_rad(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 255);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 0);
    }
    { //green
        const Color c = { 120 * (180. / M_PI), 1., 0.5 };
        const Color c_rgb = GL::Detail::HSL_to_RGB_rad(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 255);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 0);
    }
    { //blue
        const Color c = { 240. * (180. / M_PI), 1., 0.5 };
        const Color c_rgb = GL::Detail::HSL_to_RGB_rad(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 255);
    }
    {
        const Color c = { 180. * (180. / M_PI), 1., 0.1 };
        const Color c_rgb = GL::Detail::HSL_to_RGB_rad(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 51);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 51);
    }
    {
        const Color c = { 200. * (180. / M_PI), 0.5, 0.88 };
        const Color c_rgb = GL::Detail::HSL_to_RGB_rad(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 209);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 230);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 240);
    }
}

TEST_CASE( "complex to hsl conversion", "UtilFunctions" ) {

    using Color = Eigen::Matrix<float, 3, 1>;
    
    std::complex<float> cn (0.f, 1.f);
    Color col = GL::Detail::c_to_HSL(10.f, cn);
    const Color c_rgb = GL::Detail::HSL_to_RGB_deg(col);

    //REQUIRE(int(std::round(c_rgb(0) * 255)) == 0);

}

TEST_CASE( "Index flattening", "UtilFunctions") {
    using Index = Eigen::Matrix<Eigen::Index, -1, 1>;

    Index im (3);
    im(0) = im(1) = im(2) = 10;
    const Eigen::Index maxIdx = im(0) * im(1) * im(2);

    Eigen::Index k = 0;

    for(size_t z = 0; z < im(2); ++z){
        for(size_t y = 0; y < im(1); ++y){
            for(size_t x = 0; x < im(0); ++x){
                Index idx(3);
                idx(0) = z;
                idx(1) = y;
                idx(2) = x;
                const Eigen::Index i = Math::Hagedorn::Detail::index(idx, im);
                REQUIRE(i >= 0);
                REQUIRE(i <= maxIdx);     
                REQUIRE(k == i);
                k++;   
            }
        }
    }
}

TEST_CASE( "hyperbolic cutshape", "Hagedorn"){
    using Vector = Eigen::Matrix<size_t, 3, 1>;
    constexpr size_t K = 4;

    const Vector results[] = 
    {
        {0, 0, 0},
        {0, 0, 1},
        {0, 0, 2},
        {0, 0, 3},

        {0, 1, 0},
        {0, 1, 1},
        {0, 2, 0},
        {0, 3, 0},

        {1, 0, 0},
        {1, 0, 1},
        {1, 1, 0},
        {2, 0, 0},
        
        {3, 0, 0}
    };

    const auto hcs = IO::hyperbolicCutShape(3, K);

    REQUIRE(sizeof(results) / sizeof(Vector) == hcs.size());

    for(size_t i = 0; i < hcs.size(); ++i)
        for(size_t j = 0; j < 3; ++j)
            REQUIRE(results[i](j) == hcs[i](j));

}

template<class T>
[[nodiscard]] Eigen::Matrix<T, -1, 1> constexpr vec(T _v1, T _v2, T _v3) noexcept{
    Eigen::Matrix<T, -1, 1> out (3);
    out(0) = _v1;
    out(1) = _v2;
    out(2) = _v3;
    return out;
}

TEST_CASE( "Function Values", "Hagedorn" ) {
    const Eigen::Matrix<double, -1, 1> grid[] = {
        vec( -2., -2., -2. ), //0
        vec( -2., -2., 0. ), //1
        vec( -2., 0., -2. ), //2
        vec( -2., 0., -0. ), //3
        vec( 0., -2., -2. ), //4
        vec( 0., -2., 0. ), //5
        vec( 0., 0., -2. ), //6
        vec( 0., 0., 0. ) //7
    };

    const Eigen::Matrix<std::complex<double>, 4, 1> psis[] = {
        { //0
            std::complex<double>( -2.5764351288030696E-4, 2.3086258535329378E-5 ),
            std::complex<double>( -2.017463661539884E-4, -7.296617474412838E-4 ) ,
            std::complex<double>( -0.0033647777455862405, 0.001606718287491552 ),
            std::complex<double>( 0.0053542062144028065, 0.003011827401759797 )
        },
        { //1
            std::complex<double>( -0.0019037423701981498, 1.705856594319653E-4 ),
            std::complex<double>( -0.0021306839429229315, -0.0022376220120880973 ),
            std::complex<double>( -0.008716921270252179, 0.005413664935008528 ),
            std::complex<double>( 0.012503059677549187, 0.015051890401318983 )
        },
        { //2
            std::complex<double>( -0.014066859171305332, 0.0012604670072158689 ),
            std::complex<double>( -0.007974469976778167, 0.004174531091377546 ),
            std::complex<double>( 0.004668885146053407, 0.011071782004597415 ),
            std::complex<double>( 0.009269953213273835, -0.005291568321759274 )
        },
        { //3
            std::complex<double>( -0.10394081155253222, 0.00931366142716928 ),
            std::complex<double>( -0.020663344988212103, 0.030374818341262003 ),
            std::complex<double>( 0.01612322413953675, 0.02887012417400934 ),
            std::complex<double>( 0.03359942211752717, -0.0051970607454832046 )
        },
        { //4
            std::complex<double>( 6.371232637838709E-4, 0. ),
            std::complex<double>( .00268385888001607, -0.0010982025708083298 ),
            std::complex<double>( -0.007294919348626058, 0.0013302910794410344 ),
            std::complex<double>( 0.01575300569585211, -0.006315080272382563 )
        },
        { //5
            std::complex<double>( 0.004707739538032813, 0. ),
            std::complex<double>( 0.007544906279154069, -0.0091188278204724 ),
            std::complex<double>( -0.019550494103090458, 0.005847683493585132 ),
            std::complex<double>( 0.0540518151274769, -8.66617388045509E-4 )
        },
        { //6
            std::complex<double>( 0.03478575154567832, 0. ),
            std::complex<double>( -0.019493977361911066, -0.02843880805256412 ),
            std::complex<double>( 0.003191225092967026, 0.023681855459629093 ),
            std::complex<double>( 0.004652496458874799, -0.02911937169654548 )
        },
        { //7
            std::complex<double>( 0.25703386961448066, 0. ),
            std::complex<double>( -0.12490250632425211, -0.06482239926588476 ),
            std::complex<double>( 0.015900014838436553, 0.06380911520337551 ),
            std::complex<double>( 0.046969874365269504, -0.08134119466712524 )
        }
    };

    const auto file = IO::getExample();

    double error = 0.;

    for( size_t i = 0; i < 8; ++i){
        for(size_t t = 0; t < 4; ++t){
            const auto phis = Math::Hagedorn::compute(
                grid[i],
                1.,
                file.k_max,
                file.p[t],
                file.q[t],
                file.Q[t],
                file.P[t]
            );

            //calculate linear combination
            std::complex<double> res (0., 0.);
 
            for(size_t k = 0; k < file.Ks.size(); ++k){ 
                const size_t idx = Math::Hagedorn::Detail::index(file.Ks[k], file.k_max);
                res += file.c_0[t](k) * phis[idx];
            } 
            res *= std::exp(std::complex<double>(0., 1.) * file.S(0));
            const double rl2 = std::abs(res);
            error += (rl2*rl2);

            std::cout << "R:" << res << " S:" << psis[i](t) << " error: (" << std::abs(psis[i](t).real() - res.real()) << ", " << std::abs(psis[i](t).imag() - res.imag()) << ")" << std::endl;

            //REQUIRE(res == psis[i](t));
        }

        std::cout << "L2 error: " << std::sqrt(error)  << std::endl;

    }



}
