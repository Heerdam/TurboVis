
//#include <TurboDorn/voxeliser.hpp>

#include <TurboDorn/griderator.hpp>

#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

using namespace TurboDorn;

TEST_CASE("test") {

    Detail::Ray ray;

    Eigen::Vector3i num_cells;

    Griderator it(std::move(ray), num_cells);
}

/*

TEST_CASE( "rgb-hsl conversion", "UtilFunctions" ) {

    using Color = Eigen::Matrix<float, 3, 1>;

    // ---------------- DEG ----------------
    { //black
        const Color c = { 0., 0., 0. };
        const Color c_rgb = Detail::HSL_to_RGB_deg(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 0);
    }
    { //white
        const Color c = { 0., 0., 1. };
        const Color c_rgb = Detail::HSL_to_RGB_deg(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 255);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 255);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 255);
    }
    { //red
        const Color c = { 0., 1., 0.5 };
        const Color c_rgb = Detail::HSL_to_RGB_deg(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 255);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 0);
    }
    { //green
        const Color c = { 120, 1., 0.5 };
        const Color c_rgb = Detail::HSL_to_RGB_deg(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 255);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 0);
    }
    { //blue
        const Color c = { 240., 1., 0.5 };
        const Color c_rgb = Detail::HSL_to_RGB_deg(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 255);
    }

    // ---------------- RAD ----------------
    { //black
        const Color c = { 0., 0., 0. };
        const Color c_rgb = Detail::HSL_to_RGB_rad(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 0);
    }
    { //white
        const Color c = { 0., 0., 1. };
        const Color c_rgb =Detail::HSL_to_RGB_rad(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 255);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 255);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 255);
    }
    { //red
        const Color c = { 0., 1., 0.5 };
        const Color c_rgb = Detail::HSL_to_RGB_rad(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 255);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 0);
    }
    { //green
        const Color c = { 120 * ( M_PI / 180.), 1., 0.5 };
        const Color c_rgb = Detail::HSL_to_RGB_rad(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 255);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 0);
    }
    { //blue
        const Color c = { 240. * (M_PI / 180.), 1., 0.5 };
        const Color c_rgb = Detail::HSL_to_RGB_rad(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 255);
    }
    {
        const Color c = { 180. * (M_PI / 180.), 1., 0.1 };
        const Color c_rgb = Detail::HSL_to_RGB_rad(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 0);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 51);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 51);
    }
    {
        const Color c = { 200. * (M_PI / 180.), 0.5, 0.88 };
        const Color c_rgb = Detail::HSL_to_RGB_rad(c);
        REQUIRE(int(std::round(c_rgb(0) * 255)) == 209);
        REQUIRE(int(std::round(c_rgb(1) * 255)) == 230);
        REQUIRE(int(std::round(c_rgb(2) * 255)) == 240);
    }
}

TEST_CASE( "complex to hsl conversion", "UtilFunctions" ) {

    using Color = Eigen::Matrix<float, 3, 1>;
    
    std::complex<float> cn (0.f, 1.f);
    Color col = Detail::c_to_HSL(10.f, cn);
    const Color c_rgb = Detail::HSL_to_RGB_deg(col);

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
                const Eigen::Index i = Detail::index(idx, im);
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

TEST_CASE( "PSI00", "Hagedorn" ){
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
            std::complex<double>( 4.2759628436366524E-4, -0.0060297147736623 ),
            std::complex<double>( -0.0017348507365289182, -0.00279245198230175 ) ,
            std::complex<double>( -0.0021267173532104246, 6.8227618326086E-4 ),
            std::complex<double>( -1.640221599278481E-5, 6.177924623080655E-4 )
        },
        { //1
            std::complex<double>( 0.0031595329328574255, -0.044553900723141665 ),
            std::complex<double>( -0.011619113072803113, -0.0067101602532227776 ),
            std::complex<double>( -0.005610469970741988, 0.0025101920159417727 ),
            std::complex<double>( -7.516389213715678E-4, 0.0018193580853363175 )
        },
        { //2
            std::complex<double>( 0.0011623272097023578, -0.016390464100037277 ),
            std::complex<double>( -0.003068618180823126, -0.007237567042780243 ),
            std::complex<double>( -0.005448637723801008, -0.0016043287803900958 ),
            std::complex<double>( -0.002363036219306998, 0.0012327875742295601 )
        },
        { //3
            std::complex<double>( 0.008588500957804251, -0.12111005872268431 ),
            std::complex<double>( -0.02497617884667305, -0.02014020498992675 ),
            std::complex<double>( -0.015389129281395773, -0.002738280154887019 ),
            std::complex<double>( -0.008432858273859735, 9.793735686233036E-4 )
        },
        { //4
            std::complex<double>( 0.014420108039292305, 0.007877740925681847 ),
            std::complex<double>( 0.017732340728101176, -0.022463734867819145 ),
            std::complex<double>( -0.02109325449900746, -0.010039166609591596 ),
            std::complex<double>( -0.005009741986883563, 0.009532377061236885 )
        },
        { //5
            std::complex<double>( 0.1065509872549717, 0.05820906963270503 ),
            std::complex<double>( 0.02057390754367418, -0.11498019525269601 ),
            std::complex<double>( -0.060735062668945114, -0.02107210459395533 ),
            std::complex<double>( -0.025749316760644556, 0.022660742765194992 )
        },
        { //6
            std::complex<double>( 0.03919791764762447, 0.0214139200075891 ),
            std::complex<double>( 0.050182407182365524, -0.04653171658681759 ),
            std::complex<double>( -0.030021576480039475, -0.05126365458065782 ),
            std::complex<double>( -0.046428837805969805, 0.0010921629782298797 )
        },
        { //7
            std::complex<double>( 0.28963561245956093, 0.15822865623408933 ),
            std::complex<double>( 0.09102740526178693, -0.2640659112791267 ),
            std::complex<double>( -0.09763959445402326, -0.13112717703052806 ),
            std::complex<double>( -0.13937474840293745, -0.04957106660448522 )
        }
    };

    FilePathResolver::resolve();

    const auto file = IO::simulation_results_phi000();
    REQUIRE(file.has_value());
    const auto inv = Detail::computeInvariants(file.value());
    double error = 0.;
    std::cout << "Hagedorn - PSI00" << std::endl;
    for(size_t t = 0; t < 4; ++t){
        std::cout << "Time Step: " << t << std::endl;
        for( size_t i = 0; i < 8; ++i){
        
            const auto res = Detail::phi_0(
                t,
                grid[i],
                inv
            );

            const double rl2 = std::abs(psis[i](t)) - std::abs(res);
            error += (rl2*rl2);

            std::cout << "Result:" << res.real() << "+i*" << res.imag() << "\t\tSolution:" << psis[i](t).real() << " +i*" << psis[i](t).imag() << "\t\t|Solution - Result|: " << std::abs(psis[i](t).real() - res.real()) << "+i*" << std::abs(psis[i](t).imag() - res.imag()) << " [" << std::abs(std::abs(res) - std::abs(psis[i](t))) << "]" << std::endl;
            //std::cout << std::abs(std::abs(res) - std::abs(psis[i](t))) << std::endl;
            //REQUIRE(res == psis[i](t));
        }

        std::cout << "L2 error over domain: " << std::sqrt(error) << std::endl << std::endl;

    }

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

    const auto file = IO::simulation_results();
    REQUIRE(file.has_value());
    const auto inv = Detail::computeInvariants(file.value());

    double error = 0.;
    std::cout << "Hagedorn - Function Values" << std::endl;
    for(size_t t = 0; t < 4; ++t){
        for( size_t i = 0; i < 8; ++i){
        
            const std::unordered_map<Eigen::Index, std::complex<double>> phis = Detail::compute(
                t,
                grid[i],
                inv
            );

            //calculate linear combination
            std::complex<double> res (0., 0.);
 
            for(size_t k = 0; k < file.value().Ks.size(); ++k){ 
                const size_t idx = Detail::index(file.value().Ks[k], file.value().k_max);
                REQUIRE(phis.contains(idx));
                const std::complex<double>& p = (*phis.find(idx)).second;
                res += file.value().c_0[t](k) * p;
            } 
            res *= std::exp(std::complex<double>(0., 1.) * file.value().S(0));
            const double rl2 = std::abs(psis[i](t)) - std::abs(res);
            error += (rl2*rl2);

            std::cout << "Result:" << res.real() << "+i*" << res.imag() << "\t\tSolution:" << psis[i](t).real() << " +i*" << psis[i](t).imag() << "\t\t|Solution - Result|: " << std::abs(psis[i](t).real() - res.real()) << "+i*" << std::abs(psis[i](t).imag() - res.imag()) << " [" << std::abs(std::abs(res) - std::abs(psis[i](t))) << "]" << std::endl;

            //REQUIRE(res == psis[i](t));
        }

        std::cout << "L2 error over domain: " << std::sqrt(error) << std::endl << std::endl;

    }

}

using Now = std::chrono::high_resolution_clock::time_point;
//375.3 gflops

TEST_CASE( "performance simulation_results_phi000", "Hagedorn" ) {

    std::cout << "performance simulation_results_phi000" << std::endl;
    using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;

    const auto file = IO::simulation_results_phi000();
    REQUIRE(file.has_value());
    const auto inv = Detail::computeInvariants(file.value());

    Vector pos (3);
    pos.setZero();

    Vector dir(3);
    dir.setZero();
    dir(2) = 1. / 1000.;

    for(size_t k = 0; k < 10; ++k){
        double mean = 0.f;
        double last = 0.f;
        for(size_t i = 0; i < 1000; ++i){

            pos + double(i) * dir;

            const Now start = std::chrono::high_resolution_clock::now();

            const std::unordered_map<Eigen::Index, std::complex<double>> phis = Detail::Hagedorn::compute(
                0,
                pos,
                inv
            );
            const Now end = std::chrono::high_resolution_clock::now();
            const double time = std::chrono::duration<double>(end - start).count();
            mean += 1./ double(i+1) * (time - last);
            last =  time;
        }
        std::cout << "Mean for 1000 hagedorn iterations: " << mean *1000 << "ms" << std::endl;
        std::cout << "File with 1500 steps: " << (mean*1920*1080*1500)/(60*60*64) << "h" << std::endl;
    }
    std::cout << std::endl;

}

TEST_CASE( "performance simulation_results_phi121", "Hagedorn" ) {

    std::cout << "performance simulation_results_phi121" << std::endl;
    using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;

    const auto file = IO::simulation_results_phi121();
    REQUIRE(file.has_value());
    const auto inv = Detail::computeInvariants(file.value());

    Vector pos (3);
    pos.setZero();

    Vector dir(3);
    dir.setZero();
    dir(2) = 1. / 1000.;

    for(size_t k = 0; k < 10; ++k){
        double mean = 0.f;
        double last = 0.f;
        for(size_t i = 0; i < 1000; ++i){

            pos + double(i) * dir;

            const Now start = std::chrono::high_resolution_clock::now();

            const std::unordered_map<Eigen::Index, std::complex<double>> phis = Detail::Hagedorn::compute(
                0,
                pos,
                inv
            );
            const Now end = std::chrono::high_resolution_clock::now();
            const double time = std::chrono::duration<double>(end - start).count();
            mean += 1./ double(i+1) * (time - last);
            last =  time;
        }
        std::cout << "Mean for 1000 hagedorn iterations: " << mean *1000 << "ms" << std::endl;
        std::cout << "File with 1500 steps: " << (mean*1920*1080*1500)/(60*60*64) << "h" << std::endl;
    }
    std::cout << std::endl;

}

TEST_CASE( "performance simulation_results_phi412", "Hagedorn" ) {

    std::cout << "performance simulation_results_phi412" << std::endl;
    using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;

    const auto file = IO::simulation_results_phi412();
    REQUIRE(file.has_value());
    const auto inv = Detail::computeInvariants(file.value());

    Vector pos (3);
    pos.setZero();

    Vector dir(3);
    dir.setZero();
    dir(2) = 1. / 1000.;

    for(size_t k = 0; k < 10; ++k){
        double mean = 0.f;
        double last = 0.f;
        for(size_t i = 0; i < 1000; ++i){

            pos + double(i) * dir;

            const Now start = std::chrono::high_resolution_clock::now();

            const std::unordered_map<Eigen::Index, std::complex<double>> phis = Detail::compute(
                0,
                pos,
                inv
            );
            const Now end = std::chrono::high_resolution_clock::now();
            const double time = std::chrono::duration<double>(end - start).count();
            mean += 1./ double(i+1) * (time - last);
            last =  time;
        }
        std::cout << "Mean for 1000 hagedorn iterations: " << mean *1000 << "ms" << std::endl;
        std::cout << "File with 1500 steps: " << (mean*1920*1080*1500)/(60*60*64) << "h" << std::endl;
    }

    std::cout << std::endl;

}

*/
