#ifndef HDF5_HPP
#define HDF5_HPP

#include <Eigen/Eigen>

#include "filepathresolver.hpp"

#include <highfive/H5File.hpp>
#include <highfive/h5easy_bits/H5Easy_Eigen.hpp>

#include <vector>
#include <filesystem>
#include <string>   
#include <iostream>
#include <regex>
#include <exception>
#include <optional>

#include <TurboDorn/turbodorn.hpp>

namespace IO {

    namespace Detail {
        struct File {
            size_t dimensions;
            size_t timesteps;
            size_t K;
            double epsilon = 1.;  
            Eigen::VectorXcd S;       
            std::vector<Eigen::VectorXcd> c_0;
            std::vector<Eigen::VectorXcd> p, q;
            std::vector<Eigen::MatrixXcd> P, Q;

            //shapefunction
            std::vector<Eigen::VectorXi> Ks; 
            Eigen::VectorXi k_max; //max k in d direction
            std::unordered_map<Eigen::Index, bool> b_Ks;
        };//File
    }

    const Detail::File loadFromFile(const std::filesystem::path& _path, size_t _dims, size_t _K) {

        using Vector = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>;
        using Matrix = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>;

        if(!std::filesystem::exists(_path))
            throw std::runtime_error("file does not exist");

        Detail::File out;
        out.dimensions = _dims;
        out.K = _K;

        try {

            //file
            HighFive::File file(_path.string(), HighFive::File::ReadOnly);

            //timesteps
            {
                const Eigen::Matrix<double, Eigen::Dynamic, 1> steps = H5Easy::load<Eigen::Matrix<double, Eigen::Dynamic, 1>>(file, "datablock_0/wavepacket/timegrid");
                out.timesteps = steps.rows();
            }

            //S
            {
                const H5Easy::DataSet dataset = file.getDataSet("datablock_0/wavepacket/Pi/S");
                std::vector<std::vector<std::vector<std::complex<double>>>> S;
                dataset.read(S);
                out.S.resize(S.size());
                for(size_t i = 0; i < S.size(); ++i){
                    out.S(i) = S[i][0][0];
                }
            }
        
            //c_0
            {
                const Matrix tmp = H5Easy::load<Matrix>(file, "datablock_0/wavepacket/coefficients/c_0");
                for(size_t t = 0; t < out.timesteps; ++t){
                    Vector c_0(tmp.cols());
                    for(size_t d = 0; d < tmp.cols(); ++d)
                        c_0(d) = tmp(t, d);
                    out.c_0.push_back(std::move(c_0));
                }
            }
        
            //p
            {
                const H5Easy::DataSet dataset = file.getDataSet("datablock_0/wavepacket/Pi/p");
                std::vector<std::vector<std::vector<std::complex<double>>>> p;
                dataset.read(p);
                for(size_t t = 0; t < p.size(); ++t){
                    Vector temp(out.dimensions);
                    for(size_t y = 0; y < p[t].size(); ++y){
                        temp(y) = p[t][y][0];
                    }
                    out.p.push_back(std::move(temp));
                }
            }
        
            //q
            {
                const H5Easy::DataSet dataset = file.getDataSet("datablock_0/wavepacket/Pi/q");
                std::vector<std::vector<std::vector<std::complex<double>>>> q;
                dataset.read(q);
                for(size_t t = 0; t < q.size(); ++t){
                    Vector temp(out.dimensions);
                    for(size_t y = 0; y < q[t].size(); ++y){
                        temp(y) = q[t][y][0];
                    }
                    out.q.push_back(std::move(temp));
                }
            }
    
            //P
            {
                const H5Easy::DataSet dataset = file.getDataSet("datablock_0/wavepacket/Pi/P");
                std::vector<std::vector<std::vector<std::complex<double>>>> P;
                dataset.read(P);
                for(size_t t = 0; t < P.size(); ++t){
                    Matrix temp(out.dimensions, out.dimensions);
                    for(size_t x = 0; x < out.dimensions; ++x){
                        for(size_t y = 0; y < out.dimensions; ++y)
                            temp(x, y) = P[t][x][y];
                    }
                    out.P.push_back(std::move(temp));
                }
            }

            //Q
            {
                const H5Easy::DataSet dataset = file.getDataSet("datablock_0/wavepacket/Pi/Q");
                std::vector<std::vector<std::vector<std::complex<double>>>> Q;
                dataset.read(Q);
                for(size_t t = 0; t < Q.size(); ++t){
                    Matrix temp(out.dimensions, out.dimensions);
                    for(size_t x = 0; x < out.dimensions; ++x){
                        for(size_t y = 0; y < out.dimensions; ++y)
                            temp(x, y) = Q[t][x][y];
                    }
                    //std::cout << temp << std::endl;
                    out.Q.push_back(std::move(temp));
                }
            }
            

        } catch(const std::exception& _e){
            throw _e;
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
                const Eigen::Index ii = TurboDorn::Detail::index(c, out.k_max);
                //std::cout << c(0) << ", " << c(1) << ", " << c(2) << " | " << out.k_max(0) << ", " << out.k_max(1) << ", " << out.k_max(2) << " | " << ii << std::endl;
                out.b_Ks.insert( {ii, true} );
            }
        }

        return out ;
    }; //loadFromFile

    const Detail::File simulation_results() {
        try {
            FilePathResolver file;
            const auto path = file.example_files_path() /= "simulation_results.hdf5";
            return Detail::loadFromFile(path, 3, 4);
        } catch(const std::exception& _e){
            throw _e;
        }
    }//simulation_results

    const Detail::File simulation_results_phi000() {
        try {
            FilePathResolver file;
            const auto path = file.example_files_path() /= "simulation_results_phi000.hdf5";
            return Detail::loadFromFile(path, 3, 1);
        } catch(const std::exception& _e){
            throw _e;
        }
    }//simulation_results_phi000

    const Detail::File simulation_results_phi100() {
        try {
            FilePathResolver file;
            const auto path = file.example_files_path() /= "simulation_results_phi100.hdf5";
            return Detail::loadFromFile(path, 3, 2);
        } catch(const std::exception& _e){
            throw _e;
        }
    }//simulation_results_phi100

    const Detail::File simulation_results_phi121() {
        try {
            FilePathResolver file;
            const auto path = file.example_files_path() /= "simulation_results_phi121.hdf5";
            return Detail::loadFromFile(path, 3, 12);
        } catch(const std::exception& _e){
            throw _e;
        }
    }//simulation_results_phi121

    const Detail::File simulation_results_phi412() {
        try {
            FilePathResolver file;
            const auto path = file.example_files_path() /= "simulation_results_phi412.hdf5";
            return Detail::loadFromFile(path, 3, 30);
        } catch(const std::exception& _e){
            throw _e;
        }
    }//simulation_results_phi412

    std::vector<Eigen::VectorXi> hyperbolicCutShape(size_t _dim, size_t _K) {
            std::vector<Eigen::VectorXi> out;

            Eigen::VectorXi index(_dim);
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
        }; //hyperbolicCutShape

    template<class T>
    std::ostream& operator<<(std::ostream& _s, const Detail::File& _f) {
        _s << std::endl << "=====================================================" << std::endl;
        _s << "Dimensions: " << _f.dimensions << std::endl;
        _s << "Time steps: " << _f.timesteps << std::endl;
        _s << "K: " << _f.K << std::endl;
        _s << "Epsilon: " << _f.epsilon << std::endl << std::endl;
        _s << "S: " << std::endl << _f.S << std::endl;
        //_s << "c_0: " << std::endl;
        //for(size_t i = 0; i < _f.c_0.size(); ++i)
        //    _s << "[" << i << "]" << std::endl << _f.c_0[i];
        _s << "p: " << std::endl;
        for(size_t i = 0; i < _f.p.size(); ++i)
            _s << "[" << i << "]" << std::endl << _f.p[i]<< std::endl;
        _s << "q: " << std::endl;
        for(size_t i = 0; i < _f.q.size(); ++i)
            _s << "[" << i << "]" << std::endl << _f.q[i]<< std::endl;
        _s << "P: " << std::endl;
        for(size_t i = 0; i < _f.P.size(); ++i)
            _s << "[" << i << "]" << std::endl << _f.P[i]<< std::endl;
        _s << "Q: " << std::endl;    
        for(size_t i = 0; i < _f.Q.size(); ++i)
            _s << "[" << i << "]" << std::endl << _f.Q[i]<< std::endl;
        _s << std::endl << "=====================================================" << std::endl;
        return _s;
    }
}

#endif //HDF5_HqP
