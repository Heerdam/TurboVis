#ifndef HDF5_HPP
#define HDF5_HPP

#include <Eigen/Eigen>

#include "filepathresolver.hpp"
#include "math.hpp"

#include <highfive/H5File.hpp>
#include <highfive/h5easy_bits/H5Easy_Eigen.hpp>

#include <vector>
#include <filesystem>
#include <string>   
#include <iostream>
#include <regex>
#include <exception>

namespace IO {

    namespace Details {}

    template<class T>
    struct File {
        size_t dimensions;
        size_t timesteps;
        size_t K;
        T epsilon = 1.;  
        Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> S;       
        std::vector<Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1>> c_0;
        std::vector<Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1>> p, q;
        std::vector<Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>> P, Q;

        //shapefunction
        std::vector<Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1>> Ks; 
        Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1> k_max; //max k in d direction
        std::unordered_map<Eigen::Index, bool> b_Ks;
    }; //Options

    template <class T>
    [[nodiscard]] const File<T> loadFromFile(std::filesystem::path /*_path*/, size_t /*_dims*/, size_t /*_K*/);

    [[nodiscard]] File<double> getExample() noexcept;
    [[nodiscard]] File<double> getExample1() noexcept;
    [[nodiscard]] std::vector<Eigen::Matrix<Eigen::Index, -1, 1>> hyperbolicCutShape(size_t /*_dim*/, size_t /*_K*/) noexcept;

}

template <class T>
inline const IO::File<T> loadFromFile(std::filesystem::path _path, size_t _dims, size_t _K) {

    using Vector = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1>;
    using Matrix = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>;

    IO::File<T> out;
    out.dimensions = _dims;
    out.K = _K;

    try {

        //file
        HighFive::File file(_path.string(), HighFive::File::ReadOnly);

        //timesteps
        {
            const Vector steps = H5Easy::load<Vector>(file, "datablock_0/wavepacket/timegrid");
            out.timesteps = steps.rows();
        }

        //S
        {
            out.S = H5Easy::load<Vector>(file, "datablock_0/wavepacket/Pi/S");
        }

        //c_0
        {
            out.c_0 = H5Easy::load<Vector>(file, "datablock_0/wavepacket/coefficients/c_0");
        }
        
        //p
        {
            out.p = H5Easy::load<Vector>(file, "datablock_0/wavepacket/Pi/p");
        }

        //q
        {
            out.q = H5Easy::load<Vector>(file, "datablock_0/wavepacket/Pi/q");
        }

        //P
        {
            for(size_t i = 0; i < out.timesteps; ++i){
                Matrix P (out.dimensions, out.dimensions);

                for(size_t x = 0; x < out.dimensions; ++x){
                    for(size_t y = 0; y < out.dimensions; ++y){

                        std::vector<size_t> idx;
                        idx.push_back(x);
                        idx.push_back(i);
                        idx.push_back(y);         
                        const T val = H5Easy::load(file, "datablock_0/wavepacket/Pi/P", idx);
                        P(x, y) = val;

                    }
                }
                out.P.push_back(std::move(P));
            }
        }

        //Q
        {
            for(size_t i = 0; i < out.timesteps; ++i){
                Matrix Q (out.dimensions, out.dimensions);

                for(size_t x = 0; x < out.dimensions; ++x){
                    for(size_t y = 0; y < out.dimensions; ++y){

                        std::vector<size_t> idx;
                        idx.push_back(x);
                        idx.push_back(i);
                        idx.push_back(y);         
                        const T val = H5Easy::load(file, "datablock_0/wavepacket/Pi/Q", idx);
                        Q(x, y) = val;

                    }
                }
                out.Q.push_back(std::move(Q));
            }
        }
        

    } catch(const std::exception& _e){
        std::cout << _e.what() << std::endl;
    }

    // -------------- cut shape --------------
    {
        out.Ks = hyperbolicCutShape(out.dimensions, out.K);

        out.k_max = Eigen::Matrix<Eigen::Index, -1, 1>(out.dimensions);
        out.k_max.setZero();

        for(size_t i = 0; i < out.dimensions; ++i){
            for(const auto& k : out.Ks){
                out.k_max(i) = std::max(out.k_max(i), k(i));
            }
        }
 
        for(const auto& c : out.Ks){
            const Eigen::Index ii = Math::Hagedorn::Detail::index(c, out.k_max);
            //std::cout << c(0) << ", " << c(1) << ", " << c(2) << " | " << out.k_max(0) << ", " << out.k_max(1) << ", " << out.k_max(2) << " | " << ii << std::endl;
            out.b_Ks.insert( {ii, true} );
        }
    }

    return out;
}; //loadFromFile

inline IO::File<double> IO::getExample() noexcept{
    File<double> out;

    out.dimensions = 3;
    out.timesteps = 4;
    out.K = 4;
    out.epsilon = 1.;

    // -------------- P --------------
    {
        //step 1
        Eigen::Matrix<std::complex<double>, -1, -1> P(3, 3);
        P.setZero();
        P(0, 0) = std::complex(0., 1.);
        P(1, 1) = std::complex(0., 1.);
        P(2, 2) = std::complex(0., 1.);
        out.P.push_back(std::move(P));
    }
    {
        //step 2
        Eigen::Matrix<std::complex<double>, -1, -1> P(3, 3);
        P.setZero();
        P(0, 0) = std::complex(-0.45936268, 0.7602446);
        P(1, 1) = std::complex(-0.45936268, 0.7602446);
        P(2, 2) = std::complex(-0.45936268, 0.7602446);
        out.P.push_back(std::move(P));
    }
    {
        //step 3
        Eigen::Matrix<std::complex<double>, -1, -1> P(3, 3);
        P.setZero();
        P(0, 0) = std::complex(-0.698456, 0.15594369);
        P(1, 1) = std::complex(-0.698456, 0.15594369);
        P(2, 2) = std::complex(-0.698456, 0.15594369);
        out.P.push_back(std::move(P));
    }
    {
        //step 4
        Eigen::Matrix<std::complex<double>, -1, -1> P(3, 3);
        P.setZero();
        P(0, 0) = std::complex(-0.60263211, -0.52313389);
        P(1, 1) = std::complex(-0.60263211, -0.52313389);
        P(2, 2) = std::complex(-0.60263211, -0.52313389);
        out.P.push_back(std::move(P));
    }

    // -------------- Q --------------
    {
        //step 1
        Eigen::Matrix<std::complex<double>, -1, -1> Q(3, 3);
        Q.setZero();
        Q(0, 0) = std::complex(1., 0.);
        Q(1, 1) = std::complex(1., 0.);
        Q(2, 2) = std::complex(1., 0.);
        out.Q.push_back(std::move(Q));
    }
    {
        //step 2
        Eigen::Matrix<std::complex<double>, -1, -1> Q(3, 3);
        Q.setZero();
        Q(0, 0) = std::complex(0.7602446, 0.91872537);
        Q(1, 1) = std::complex(0.7602446, 0.91872537);
        Q(2, 2) = std::complex(0.7602446, 0.91872537);
        out.Q.push_back(std::move(Q));
    }
    {
        //step 3
        Eigen::Matrix<std::complex<double>, -1, -1> Q(3, 3);
        Q.setZero();
        Q(0, 0) = std::complex(0.15594369, 1.396912);
        Q(1, 1) = std::complex(0.15594369, 1.396912);
        Q(2, 2) = std::complex(0.15594369, 1.396912);
        out.Q.push_back(std::move(Q));
    }
    {
        //step 4
        Eigen::Matrix<std::complex<double>, -1, -1> Q(3, 3);
        Q.setZero();
        Q(0, 0) = std::complex(-0.52313389, 1.20526423);
        Q(1, 1) = std::complex(-0.52313389, 1.20526423);
        Q(2, 2) = std::complex(-0.52313389, 1.20526423);
        out.Q.push_back(std::move(Q));
    }

    // -------------- p --------------
    {
        //step 0
        Eigen::Matrix<std::complex<double>, -1, 1> p(3);
        p(0) = std::complex(1., 0.);
        p(1) = std::complex(0., 0.);
        p(2) = std::complex(0., 0.);
        out.p.push_back(std::move(p));
    }
    {
        //step 1
        Eigen::Matrix<std::complex<double>, -1, 1> p(3);
        p(0) = std::complex(0.9899259395420181, 0.);
        p(1) = std::complex(0.2296813424663948, 0.);
        p(2) = std::complex(0., 0.);
        out.p.push_back(std::move(p));
    }
    {
        //step 2
        Eigen::Matrix<std::complex<double>, -1, 1> p(3);
        p(0) = std::complex(0.5051716940837586, 0.);
        p(1) = std::complex(0.34922799931828014, 0.);
        p(2) = std::complex(0., 0.);
        out.p.push_back(std::move(p));
    }
    {
        //step 3
        Eigen::Matrix<std::complex<double>, -1, 1> p(3);
        p(0) = std::complex(-0.22181783749679396, 0.);
        p(1) = std::complex(0.3013160567921123, 0.);
        p(2) = std::complex(0., 0.);
        out.p.push_back(std::move(p));
    }

    // -------------- q --------------
    {
        //step 0
        Eigen::Matrix<std::complex<double>, -1, 1> q(3);
        q(0) = std::complex(-0.5, 0.);
        q(1) = std::complex(-0.5, 0.);
        q(2) = std::complex(0., 0.);
        out.q.push_back(std::move(q));
    }
    {
        //step 1
        Eigen::Matrix<std::complex<double>, -1, 1> q(3);
        q(0) = std::complex(0.5386030713277431, 0. );
        q(1) = std::complex(-0.3801222985378165, 0.);
        q(2) = std::complex(0., 0.);
        out.q.push_back(std::move(q));
    }
    {
        //step 2
        Eigen::Matrix<std::complex<double>, -1, 1> q(3);
        q(0) = std::complex(1.3189401498906073, 0.);
        q(1) = std::complex(-0.07797184738268172, 0.);
        q(2) = std::complex(0., 0.);
        out.q.push_back(std::move(q));
    }
    {
        //step 3
        Eigen::Matrix<std::complex<double>, -1, 1> q(3);
        q(0) = std::complex(1.4668311743131408, 0.);
        q(1) = std::complex(0.26156694714445683, 0.);
        q(2) = std::complex(0., 0.);
        out.q.push_back(std::move(q));
    }

    // -------------- c_0 --------------
    {
        //step 1
        Eigen::Matrix<std::complex<double>, -1, 1> c_0(13);
        c_0.setZero();
        c_0(0) = std::complex<double>(1., 0.);
        out.c_0.push_back(std::move(c_0));
    }
    {
        //step 2
        Eigen::Matrix<std::complex<double>, -1, 1> c_0(13);
        c_0(0) = std::complex<double>(1.0, 2.33955611e-18);
        c_0(1) = std::complex<double>(8.42028803e-35, 2.18267605e-35);
        c_0(2) = std::complex<double>(-1.17869486e-19, 3.15323351e-19);
        c_0(3) = std::complex<double>(-1.37868073e-36, 1.76751386e-35);
        c_0(4) = std::complex<double>(-1.64234376e-18, 1.56297375e-18);
        c_0(5) = std::complex<double>(3.06694315e-35, 1.88248259e-36);
        c_0(6) = std::complex<double>(-2.00649974e-18, -3.68142805e-19);
        c_0(7) = std::complex<double>(-2.64217875e-18, -1.19804980e-18);
        c_0(8) = std::complex<double>(-1.96300292e-18, 1.00000000e+00);
        c_0(9) = std::complex<double>(-2.78458174e-35, 1.53115593e-35);
        c_0(10) = std::complex<double>(-1.90193417e-18, -2.37678681e-18);
        c_0(11) = std::complex<double>(2.00602071e-18, -1.99652006e-20);
        c_0(12) = std::complex<double>(1.83157351e-18, -3.32255209e-18);
        out.c_0.push_back(std::move(c_0));
    }
    {
        //step 3
        Eigen::Matrix<std::complex<double>, -1, 1> c_0(13);
        c_0(0) = std::complex<double>(1.00000000e+00, -5.11438227e-18);
        c_0(1) = std::complex<double>(1.28108837e-34, 2.91799816e-35);
        c_0(2) = std::complex<double>(4.87128610e-19, 6.74850601e-19);
        c_0(3) = std::complex<double>(-2.53780683e-35, -8.11616688e-36);
        c_0(4) = std::complex<double>(-1.81928745e-18, 3.16583998e-18);
        c_0(5) = std::complex<double>(2.45069175e-35, -9.18102116e-35);
        c_0(6) = std::complex<double>(-3.27780635e-18, -4.47067340e-19);
        c_0(7) = std::complex<double>(-4.54664503e-18, -2.21372383e-18);
        c_0(8) = std::complex<double>(1.18928369e-17, 1.00000000e+00);
        c_0(9) = std::complex<double>(-6.06717711e-35, 1.28989687e-34);
        c_0(10) = std::complex<double>(-3.41241433e-18, -4.68812046e-18);
        c_0(11) = std::complex<double>(4.08075704e-18, 6.68898666e-18);
        c_0(12) = std::complex<double>(2.27839682e-18, -5.76631175e-18);
        out.c_0.push_back(std::move(c_0));
    }
    {
        //step 4
        Eigen::Matrix<std::complex<double>, -1, 1> c_0(13);
        c_0(0) = std::complex<double>(1.00000000e+00, -8.91882937e-19);
        c_0(1) = std::complex<double>(2.86451329e-34, 3.20344182e-35);
        c_0(2) = std::complex<double>(-3.03652121e-19, -7.92773766e-19);
        c_0(3) = std::complex<double>(4.48674371e-35, 3.76248537e-36);
        c_0(4) = std::complex<double>(-5.08593138e-18, 2.05724300e-18);
        c_0(5) = std::complex<double>(-6.36527376e-36, -2.47841399e-35);
        c_0(6) = std::complex<double>(-3.36294256e-18, 2.09736066e-19);
        c_0(7) = std::complex<double>(-8.33621159e-18, -2.18035017e-18);
        c_0(8) = std::complex<double>(2.34408811e-19, 1.00000000e+00);
        c_0(9) = std::complex<double>(-3.56632043e-35, 2.81898695e-34);
        c_0(10) = std::complex<double>(-9.44547227e-19, -1.03208891e-17);
        c_0(11) = std::complex<double>(5.18094648e-18, -1.03681495e-17);
        c_0(12) = std::complex<double>(1.06469989e-17, -6.54653106e-18);
        out.c_0.push_back(std::move(c_0));
    }

    // -------------- S --------------
    {
        out.S.resize(4);
        out.S(0) = { 0., 0. };
        out.S(1) = { 0.47293507579737626, 0. };
        out.S(2) = { 0.5695306388253558, 0. };
        out.S(3) = { 0.12672250102084176, 0. };
    }

    // -------------- cut shape --------------
    {
        out.Ks = hyperbolicCutShape(out.dimensions, out.K);

        out.k_max = Eigen::Matrix<Eigen::Index, -1, 1>(out.dimensions);
        out.k_max.setZero();

        for(size_t i = 0; i < out.dimensions; ++i){
            for(const auto& k : out.Ks){
                out.k_max(i) = std::max(out.k_max(i), k(i));
            }
        }
 
        for(const auto& c : out.Ks){
            const Eigen::Index ii = Math::Hagedorn::Detail::index(c, out.k_max);
            //std::cout << c(0) << ", " << c(1) << ", " << c(2) << " | " << out.k_max(0) << ", " << out.k_max(1) << ", " << out.k_max(2) << " | " << ii << std::endl;
            out.b_Ks.insert( {ii, true} );
        }
    }

    return out;
}; //IO::getExample

inline IO::File<double> IO::getExample1() noexcept{
       File<double> out;

    out.dimensions = 3;
    out.timesteps = 4;
    out.K = 30;
    out.epsilon = 1.;

    // -------------- P --------------
    {
        //step 1
        Eigen::Matrix<std::complex<double>, -1, -1> P(3, 3);
        P.setZero();
        P(0, 0) = std::complex(0., 1.);
        P(1, 1) = std::complex(0., 1.);
        P(2, 2) = std::complex(0., 1.);
        out.P.push_back(std::move(P));
    }
    {
        //step 2
        Eigen::Matrix<std::complex<double>, -1, -1> P(3, 3);
        P.setZero();
        P(0, 0) = std::complex(-0.45936268, 0.7602446);
        P(1, 1) = std::complex(-0.45936268, 0.7602446);
        P(2, 2) = std::complex(-0.45936268, 0.7602446);
        out.P.push_back(std::move(P));
    }
    {
        //step 3
        Eigen::Matrix<std::complex<double>, -1, -1> P(3, 3);
        P.setZero();
        P(0, 0) = std::complex(-0.698456, 0.15594369);
        P(1, 1) = std::complex(-0.698456, 0.15594369);
        P(2, 2) = std::complex(-0.698456, 0.15594369);
        out.P.push_back(std::move(P));
    }
    {
        //step 4
        Eigen::Matrix<std::complex<double>, -1, -1> P(3, 3);
        P.setZero();
        P(0, 0) = std::complex(-0.60263211, -0.52313389);
        P(1, 1) = std::complex(-0.60263211, -0.52313389);
        P(2, 2) = std::complex(-0.60263211, -0.52313389);
        out.P.push_back(std::move(P));
    }

    // -------------- Q --------------
    {
        //step 1
        Eigen::Matrix<std::complex<double>, -1, -1> Q(3, 3);
        Q.setZero();
        Q(0, 0) = std::complex(1., 0.);
        Q(1, 1) = std::complex(1., 0.);
        Q(2, 2) = std::complex(1., 0.);
        out.Q.push_back(std::move(Q));
    }
    {
        //step 2
        Eigen::Matrix<std::complex<double>, -1, -1> Q(3, 3);
        Q.setZero();
        Q(0, 0) = std::complex(0.7602446, 0.91872537);
        Q(1, 1) = std::complex(0.7602446, 0.91872537);
        Q(2, 2) = std::complex(0.7602446, 0.91872537);
        out.Q.push_back(std::move(Q));
    }
    {
        //step 3
        Eigen::Matrix<std::complex<double>, -1, -1> Q(3, 3);
        Q.setZero();
        Q(0, 0) = std::complex(0.15594369, 1.396912);
        Q(1, 1) = std::complex(0.15594369, 1.396912);
        Q(2, 2) = std::complex(0.15594369, 1.396912);
        out.Q.push_back(std::move(Q));
    }
    {
        //step 4
        Eigen::Matrix<std::complex<double>, -1, -1> Q(3, 3);
        Q.setZero();
        Q(0, 0) = std::complex(-0.52313389, 1.20526423);
        Q(1, 1) = std::complex(-0.52313389, 1.20526423);
        Q(2, 2) = std::complex(-0.52313389, 1.20526423);
        out.Q.push_back(std::move(Q));
    }

    // -------------- p --------------
    {
        //step 0
        Eigen::Matrix<std::complex<double>, -1, 1> p(3);
        p(0) = std::complex(1., 0.);
        p(1) = std::complex(0., 0.);
        p(2) = std::complex(0., 0.);
        out.p.push_back(std::move(p));
    }
    {
        //step 1
        Eigen::Matrix<std::complex<double>, -1, 1> p(3);
        p(0) = std::complex(0.9899259395420181, 0.);
        p(1) = std::complex(0.2296813424663948, 0.);
        p(2) = std::complex(0., 0.);
        out.p.push_back(std::move(p));
    }
    {
        //step 2
        Eigen::Matrix<std::complex<double>, -1, 1> p(3);
        p(0) = std::complex(0.5051716940837586, 0.);
        p(1) = std::complex(0.34922799931828014, 0.);
        p(2) = std::complex(0., 0.);
        out.p.push_back(std::move(p));
    }
    {
        //step 3
        Eigen::Matrix<std::complex<double>, -1, 1> p(3);
        p(0) = std::complex(-0.22181783749679396, 0.);
        p(1) = std::complex(0.3013160567921123, 0.);
        p(2) = std::complex(0., 0.);
        out.p.push_back(std::move(p));
    }

    // -------------- q --------------
    {
        //step 0
        Eigen::Matrix<std::complex<double>, -1, 1> q(3);
        q(0) = std::complex(-0.5, 0.);
        q(1) = std::complex(-0.5, 0.);
        q(2) = std::complex(0., 0.);
        out.q.push_back(std::move(q));
    }
    {
        //step 1
        Eigen::Matrix<std::complex<double>, -1, 1> q(3);
        q(0) = std::complex(0.5386030713277431, 0. );
        q(1) = std::complex(-0.3801222985378165, 0.);
        q(2) = std::complex(0., 0.);
        out.q.push_back(std::move(q));
    }
    {
        //step 2
        Eigen::Matrix<std::complex<double>, -1, 1> q(3);
        q(0) = std::complex(1.3189401498906073, 0.);
        q(1) = std::complex(-0.07797184738268172, 0.);
        q(2) = std::complex(0., 0.);
        out.q.push_back(std::move(q));
    }
    {
        //step 3
        Eigen::Matrix<std::complex<double>, -1, 1> q(3);
        q(0) = std::complex(1.4668311743131408, 0.);
        q(1) = std::complex(0.26156694714445683, 0.);
        q(2) = std::complex(0., 0.);
        out.q.push_back(std::move(q));
    }

    // -------------- c_0 --------------
    {
        const auto splitter = [](const std::string& _in, std::vector<std::string>& _out, const std::string& _delimeter) {
            const std::regex r(_delimeter);
            std::sregex_token_iterator begin(_in.begin(), _in.end(), r, -1);
            if (std::distance(begin, std::sregex_token_iterator()) == 0) return;
            for (auto it = begin; it != std::sregex_token_iterator(); ++it) {
                if((*it).str().empty()) continue;
                _out.push_back(*it);
            }
        };

        const auto path = FilePathResolver::ASSETDIR() + "c0.txt";
        std::ifstream file (path);
        std::string line;
        for(size_t i = 0; i < 4; ++i) {

            std::getline(file, line);
            std::vector<std::string> vals;
            splitter(line, vals, "\t");

            Eigen::Matrix<std::complex<double>, -1, 1> c_0(vals.size());
            for(size_t k = 0; k < vals.size(); ++k){
                c_0(k) = std::stod(vals[k]);
            }
            out.c_0.push_back(std::move(c_0));
        }
 
    }

    // -------------- S --------------
    {
        out.S.resize(4);
        out.S(0) = { 0, 0 };
        out.S(1) = { 0.47293507579737626, 0 };
        out.S(2) = { 0.5695306388253558, 0 };
        out.S(3) = { 0.12672250102084176, 0 };
    }

    // -------------- cut shape --------------
    {
        out.Ks = hyperbolicCutShape(out.dimensions, out.K);

        out.k_max = Eigen::Matrix<Eigen::Index, -1, 1>(out.dimensions);
        out.k_max.setZero();

        for(size_t i = 0; i < out.dimensions; ++i){
            for(const auto& k : out.Ks){
                out.k_max(i) = std::max(out.k_max(i), k(i));
            }
        }
 
        for(const auto& c : out.Ks){
            const Eigen::Index ii = Math::Hagedorn::Detail::index(c, out.k_max);
            //std::cout << c(0) << ", " << c(1) << ", " << c(2) << " | " << out.k_max(0) << ", " << out.k_max(1) << ", " << out.k_max(2) << " | " << ii << std::endl;
            out.b_Ks.insert( {ii, true} );
        }
    }

    return out;
} //IO::getExample1

inline std::vector<Eigen::Matrix<Eigen::Index, -1, 1>> IO::hyperbolicCutShape(size_t _dim, size_t _K) noexcept {
    std::vector<Eigen::Matrix<Eigen::Index, -1, 1>> out;

    Eigen::Matrix<Eigen::Index, -1, 1> index(_dim);
    index.setZero();

    while (true) {
        for (index(_dim - 1) = 0; index(_dim - 1) <= _K; ++index(_dim - 1)) {
            Eigen::Index p = 1;
            for (Eigen::Index d = 0; d < _dim; ++d)
                p *= (1 + index(d));
            if (p <= _K)
                out.push_back(index);
        }

        bool done = false;
        for (Eigen::Index d = _dim - 2; d >= 0; --d) {
            index(d) += 1;

            if (index(d) >= _K) {
                if (d == 0)
                    done = true;
                else
                    index(d) = 0;
            } else
                break;
        }
        if (done) break;
    }

    return out;
}; //IO::hyperbolicCutShape

#endif //HDF5_HqP