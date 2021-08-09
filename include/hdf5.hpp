#ifndef HDF5_HPP
#define HDF5_HPP

#include <Eigen/Eigen>

#include <highfive/H5File.hpp>
#include <highfive/h5easy_bits/H5Easy_Eigen.hpp>

#include <vector>
#include <filesystem>
#include <string>
#include <iostream>

namespace IO {

    namespace Details {

        template<class Node>
        std::vector<unsigned char> extractOpaque(const Node& _node){
            char* data = new char[_node.getStorageSize()];
            _node.read(data, _node.getDataType());
            std::vector<unsigned char> out;
            for (size_t i = 0; i < _node.getStorageSize(); ++i){
                if(data[i] < 0x41 || data[i] > 0x7A) continue;
                //std::cout << size_t(data[i]) << " ";
                out.push_back(data[i]);
            }


            //std::cout << std::endl;
            std::string ss (out.begin(), out.end());
            std::cout << ss << std::endl;
            return out;
        }

    }

    template<class T>
    struct File {
        size_t dimensions;
        size_t K;
        Eigen::Matrix<Eigen::Index, -1, 1> k_max;
        Eigen::Matrix<std::complex<T>, -1, 1> S;
        std::vector<Eigen::Matrix<Eigen::Index, -1, 1>> Ks;
        std::vector<Eigen::Matrix<std::complex<T>, -1, 1>> c_0;
        std::vector<Eigen::Matrix<std::complex<T>, -1, 1>> p, q;
        std::vector<Eigen::Matrix<std::complex<T>, -1, -1>> P, Q;
    }; //Options

    template <class T, size_t Dim>
    [[nodiscard]] const File<T> loadFromFile(std::filesystem::path _path) {
        File<T> out;
        using Matrix = Eigen::Matrix<std::complex<double>, -1, -1>;
        //file
        HighFive::File file(_path.string(), HighFive::File::ReadOnly);

        //basishapes
        {
            const Eigen::Matrix<Eigen::Index, -1, 1> hashes = H5Easy::load<Eigen::Matrix<Eigen::Index, -1, 1>>(file, "datablock_0/wavepacket/basis_shape_hash");
            std::cout << hashes << std::endl
                      << std::endl;
            for (Eigen::Index i = 0; i < hashes.rows(); ++i) {
                std::stringstream ss;
                ss << "/datablock_0/wavepacket/basisshapes/basis_shape_" << hashes(i);
                const auto dataset = file.getDataSet(ss.str());
                auto dd = dataset.listAttributeNames();
                const auto a_k = dataset.getAttribute("K");
                const auto a_d = dataset.getAttribute("dimension");
                const auto a_t = dataset.getAttribute("type");
                

                const auto test1 = Details::extractOpaque(a_k);
                const auto test2 = Details::extractOpaque(a_d);
                const auto test3 = Details::extractOpaque(a_t);

                auto sf = a_k.getDataType();
                
            }
        }

        Matrix c_0 = H5Easy::load<Matrix>(file, "datablock_0/wavepacket/coefficients/c_0");

        return out;
    }; //loadFromFile

    [[nodiscard]] File<double> getExample() noexcept;
    [[nodiscard]] std::vector<Eigen::Matrix<Eigen::Index, -1, 1>> hyperbolicCutShape(size_t /*_dim*/, size_t /*_K*/) noexcept;

}

inline IO::File<double> IO::getExample() noexcept{
    File<double> out;

    out.dimensions = 3;
    out.K = 4;

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
        //step 1
        Eigen::Matrix<std::complex<double>, -1, 1> p(3);
        p(0) = std::complex(1., 0.);
        p(1) = std::complex(0., 0.);
        p(2) = std::complex(0., 0.);
        out.p.push_back(std::move(p));
    }
    {
        //step 1
        Eigen::Matrix<std::complex<double>, -1, 1> p(3);
        p(0) = std::complex(0.7602446, 0.);
        p(1) = std::complex(-0.45936268, 0.);
        p(2) = std::complex(0., 0.);
        out.p.push_back(std::move(p));
    }
    {
        //step 1
        Eigen::Matrix<std::complex<double>, -1, 1> p(3);
        p(0) = std::complex(0.15594369, 0.);
        p(1) = std::complex(-0.698456, 0.);
        p(2) = std::complex(0., 0.);
        out.p.push_back(std::move(p));
    }
    {
        //step 1
        Eigen::Matrix<std::complex<double>, -1, 1> p(3);
        p(0) = std::complex(-0.52313389, 0.);
        p(1) = std::complex(-0.60263211, 0.);
        p(2) = std::complex(0., 0.);
        out.p.push_back(std::move(p));
    }

    // -------------- q --------------
    {
        //step 1
        Eigen::Matrix<std::complex<double>, -1, 1> q(3);
        q(0) = std::complex(0., 0.);
        q(1) = std::complex(1., 0.);
        q(2) = std::complex(0., 0.);
        out.q.push_back(std::move(q));
    }
    {
        //step 1
        Eigen::Matrix<std::complex<double>, -1, 1> q(3);
        q(0) = std::complex(0.9187253698655619, 0.);
        q(1) = std::complex(0.760244597075633, 0.);
        q(2) = std::complex(0., 0.);
        out.q.push_back(std::move(q));
    }
    {
        //step 1
        Eigen::Matrix<std::complex<double>, -1, 1> q(3);
        q(0) = std::complex(1.3969119972732962, 0.);
        q(1) = std::complex(0.15594369476536343, 0.);
        q(2) = std::complex(0., 0.);
        out.q.push_back(std::move(q));
    }
    {
        //step 1
        Eigen::Matrix<std::complex<double>, -1, 1> q(3);
        q(0) = std::complex(1.205264227168689, 0.);
        q(1) = std::complex(-0.5231338942889137, 0.);
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
        out.S(1) = { 0.17461399965914637, 0. };
        out.S(2) = { 0.05445990452930107, 0. };
        out.S(3) = { -0.15762864220149592, 0. };
    }

    out.Ks = hyperbolicCutShape(out.dimensions, out.K);

    out.k_max = Eigen::Matrix<Eigen::Index, -1, 1>(out.dimensions);
    out.k_max.setZero();

    for(size_t i = 0; i < out.dimensions; ++i){
        for(const auto& k : out.Ks){
            out.k_max(i) = std::max(out.k_max(i), k(i));
        }
    }

    return out;
}; //IO::getExample

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