#ifndef TURBODORN_HPP
#define TURBODORN_HPP

#include <complex>
#include <numeric>
#include <chrono>
#include <vector>
#include <filesystem>
#include <fstream>
#include <cassert>
#include <exception>
#include <optional>
#include <unordered_map>

#include <Eigen/Eigen>

#include "griderator.hpp"

namespace TurboDorn {

    namespace Detail {

        struct Invariants {
            size_t dimensions;
            Eigen::VectorXi k; //k extends
            std::unordered_map<Eigen::Index, bool> k_shape; //lookup for shape
            std::vector<Eigen::VectorXcd> p;
            std::vector<Eigen::VectorXcd> q;
            //phi_0
            std::vector<std::complex<double>> pre;
            std::complex<double> i_2_E_2;
            std::vector<Eigen::MatrixXcd> P_Q_1;
            std::vector<Eigen::RowVectorXcd> i_E_2_p;
            //phi
            //sqrt*Q-1
            std::vector<Eigen::MatrixXcd> Q_1;
            //Q-1*QT
            std::vector<Eigen::MatrixXcd> Q_1_Q_T;

            Invariants() = default;
            Invariants(Invariants&& _in) = default;               
            [[nodiscard]] Invariants& operator=(Invariants&& _in) noexcept = default;
        };//Invariants

        std::unique_ptr<Detail::Invariants> computeInvariants(const File& _file) {
            std::unique_ptr<Detail::Invariants> out = std::make_unique();
            out->dimensions = _file.dimensions;
            out->k = _file.k_max;
            out->i_2_E_2 = std::complex<double>(0., 1.) / (2. * _file.epsilon * _file.epsilon);
            out->p = _file.p;
            out->q = _file.q;
            out->k_shape = _file.b_Ks;

            for(size_t t = 0; t < _file.timesteps; ++t){
                //phi0
                out->pre.push_back( std::pow(double(M_PI) * _file.epsilon * _file.epsilon, - double(_file.dimensions) / 4.) * std::pow(_file.Q[t].determinant(), -0.5) );
                out->P_Q_1.push_back( _file.P[t] * _file.Q[t].inverse() );

                out->i_E_2_p.push_back( (std::complex<double>(0., 1.) / _file.epsilon * _file.epsilon) * _file.p[t].transpose() );
                //phi
                out->Q_1.push_back( std::sqrt(2. / (_file.epsilon * _file.epsilon)) * _file.Q[t].inverse() );
                out->Q_1_Q_T.push_back( _file.Q[t].inverse() * _file.Q[t].conjugate() );
            }
            return out;
        }

        /*
        i: index
        e: extents, # units
        */
        inline Eigen::Index index(const Eigen::VectorXi& _i, const Eigen::VectorXi& _e ) {
            assert(_i.size() == _e.size());
            Eigen::Index out = _i(0);
            for (Eigen::Index k = 1; k < _i.size(); ++k) {
                out *= _e(k);
                out += _i(k);
            }
            return out;
        }; //index

        struct StreamingOctreeOptions {

        };

        /*
            voxel size: size of a voxel in world coords
            depth: how many subdivides the octree has
            chunk size: symmetric 2^n extends -> size: 4*2*(2^n)^3 bytes (n = 32 (262'144 bytes) or 64 (2'097'152 bytes))

            only lowest level has data. rest only boolmaps

            format:
            double                              voxel size
            size_t                              depth
            size_t                              chunk size
            n bits with 1 byte/subdivision      ((8^(n+1))-1)/7-1 nodes
            chunks...                           chunk_size^3 bytes for every none-empty chunk

        */

        template<class T = double>
        class Chunkerator {

            

            iVec extends;

            bVec bmap;
            cVec temp_chunk;

            //iterator state
            const bool isEndIterator = false;

            //file meta information
            std::filesystem::path path;
            dVec camPos, camPlaneNorm, halfExtCamPlane;


        public:

            using iterator_category = std::forward_iterator_tag;

            /*
                path - path to the file
                p - position of the camera plane in world coords
                n - normalized plane normal of the camera plane
                e - half extends of the camera plane in world units
            */
            Chunkerator(std::filesystem::path&& _path, dVec&& _p, dVec&& _n, dVec&& _e) : 
                path(std::move(_path)), camPos(std::move(_p)),  camPlaneNorm(std::move(_n)), halfExtCamPlane(std::move(_e)){}

            //creates an end chunkerator
            Chunkerator() : isEndIterator(true) {};

            void intitialise() {
                if(!isEndIterator) throw std::runtime_error("Can't initialize an end chunkerator!");
                
            }

            void operator*() {

            }

            void operator++() {

            }

            void operator++(int) {

            }

            bool operator==(const Chunkerator<T>& _other) const noexcept {
                if(isEndIterator && _other.isEndIterator) return true;
            }

        };//Chunkerator

















        template<class Vec, class T>
        [[nodiscard]] bool intersect(const Vec& _r_o, const Vec& _r_d, const Vec& _low, const Vec& _high, T _tmax, T& _t) noexcept {
            _t = -std::numeric_limits<T>::infinity();
            for (size_t i = 0; i < _r_o.size(); ++i) {
                if (std::abs(_r_d[i]) < std::numeric_limits<T>::epsilon()){
                    if (_r_o[i] < _low[i] || _r_o[i] > _high[i]) 
                        return false;
                } else {
                    const T ood = 1. / _r_d[i];
                    T t1 = (_low[i] - _r_o[i]) * ood;
                    T t2 = (_high[i] - _r_o[i]) * ood;
                    if (t1 > t2) std::swap(t1, t2);
                    _t = std::max(_t, t1);
                    _tmax = std::min(_tmax, t2);
                    if (_t > _tmax) return false;
                }
            }
            return true;
        }//intersect

        template <class T>
        [[nodiscard]] Eigen::Matrix<T, 3, 1> c_to_HSL(T _max, const std::complex<T>& _c) noexcept {
            const T H = std::clamp(std::abs(std::fmod(std::arg(_c), 2. * M_PI)), 0., 2. * M_PI);
            const T S = 1.;
            const T L = std::clamp(std::abs(_max * std::atan(std::abs(_c)) / (0.5 * M_PI)), 0., 1.);
            return { H, S, L };
        }//c_to_HSL

        /*
            h: [0, 2pi]
            s: [0, 1]
            l: [0, 1]
            rgb: [0, 1]
        */
        template <class T>
        [[nodiscard]] Eigen::Matrix<T, 3, 1> HSL_to_RGB_rad(const Eigen::Matrix<T, 3, 1>& _hsl) noexcept {
            return HSL_to_RGB_deg<T>( { _hsl(0) * T( 180. / M_PI), _hsl(1), _hsl(2) } );
        }; //HSL_to_RGB_rad

        /*.
            h: [0, 360]
            s: [0, 1]
            l: [0, 1]
            rgb: [0, 1]
        */
        template <class T>
        [[nodiscard]] Eigen::Matrix<T, 3, 1> HSL_to_RGB_deg(const Eigen::Matrix<T, 3, 1>& _hsl) noexcept {

            const T H = _hsl(0);
            const T S = _hsl(1);
            const T L = _hsl(2);

            assert(0. <= H && H <= 360.);
            assert(0. <= S && S <= 1.);
            assert(0. <= L && L <= 1.);

            const T C = ( T(1.) - std::abs( T(2.) * L - T(1.) ) ) * S;
            const T X = C * (T(1.) - std::abs(std::fmod(H / T(60.), T(2.)) - T(1.)));
            const T m = L - C * T(0.5);

            switch(size_t(H / 60.)){
                case 0: return { C + m, X + m, m};
                case 1: return { X + m, C + m, m};
                case 2: return { m, C + m, X + m};
                case 3: return { m, X + m, C + m};
                case 4: return { X + m, m, C + m};
                case 5: return { C + m, m, X + m};
                default: return { 0., 0., 0.};
            }

        }; //HSL_to_RGB_deg

        template <class T>
        [[nodiscard]] Eigen::Matrix<T, 3, 1> rgb_to_gs(const Eigen::Matrix<T, 3, 1>& _rgb) noexcept {
            const T gs = (_rgb(0) + _rgb(1) + _rgb(2)) / 3.;
            return { gs, gs, gs };
        }; //rgb_to_gs

        template <class T>
        inline std::complex<T> phi_0 (
            size_t _t, 
            const Eigen::Matrix<T, Eigen::Dynamic, 1>& _x, 
            const Detail::Invariants<T>& _inv
        ) noexcept {
            const Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> xq = _x - _inv.q[_t];
            const Eigen::Matrix<std::complex<T>, 1, Eigen::Dynamic> xqt = xq.transpose();
            const std::complex<T> e1 = _inv.i_2_E_2 * xqt * _inv.P_Q_1[_t] * xq;
            const std::complex<T> e2 = _inv.i_E_2_p[_t] * xq;
            return _inv.pre[_t] * std::exp(e1 + e2);
        }; //phi_0

        template <class T>
        inline Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> phi (
            size_t _t,
            const Eigen::Matrix<T, Eigen::Dynamic, 1>& _x,
            const std::unordered_map<Eigen::Index, std::complex<T>>& _phis,
            const Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1>& _index,
            const Detail::Invariants<T>& _inv
        ) noexcept {

            using Index = Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1>;
            using Vector = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1>;
            using Matrix = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>;

            const auto xq = _x - _inv.q[_t];

            Vector res (_inv.dimensions);

            Vector kp (_inv.dimensions);
            for (size_t j = 0; j < _inv.dimensions; ++j) {
                if (_index(j) - 1 < 0) {
                    kp(j) = { 0., 0. };
                    continue;
                }

                Index k_1 = _index;
                k_1(j) = std::max<Eigen::Index>(k_1(j) - 1, 0ll);
                const Eigen::Index ii = Detail::index(k_1, _inv.k);
                assert(_phis.contains(ii));
                const std::complex<double>& p = (*_phis.find(ii)).second;
                kp(j) = std::sqrt(_index(j)) * p;  
            }

            const Eigen::Index ii = Detail::index(_index, _inv.k);
            assert(_phis.contains(ii));
            const std::complex<double>& p = (*_phis.find(ii)).second;
            auto phi_t = _inv.Q_1[_t] * xq * p - _inv.Q_1_Q_T[_t] * kp;

            Vector phi (_inv.dimensions);
            for(size_t i = 0; i < _inv.dimensions; ++i){
                const T sk = std::sqrt(T(_index(i)) + 1.);
                const std::complex<T> skc = std::complex<T>(sk, 0.);
                phi(i) = phi_t(i) / skc;
            }

            return phi;
        }; //phi

        template <class T>
        inline std::unordered_map<Eigen::Index, std::complex<T>> computeCube (
            size_t _t, 
            const Eigen::Matrix<T, Eigen::Dynamic, 1>& _x, 
            const Detail::Invariants<T>& _inv
        ) noexcept {

            using Index = Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1>;
            using Vector = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1>;
            using Matrix = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>;

            const size_t dim = _inv.dimensions;

            size_t size = _inv.k(0) + 1;
            for (size_t i = 1; i < _inv.k.size(); ++i)
                size *= (_inv.k(i)+1);

            std::unordered_map<Eigen::Index, std::complex<T>> phis;

            //iterate over ks
            Index index(dim);
            index.fill(0);

            bool first = true;

            while (true) {
                for (index(dim-1) = 0; index(dim-1) <= _inv.k(dim-1); ++index(dim-1)) { 
                    if (first) {
                        first = false;
                        const auto phi0 = Detail::phi_0(_t, _x, _inv);
                        phis.insert( {0, phi0} );
                        --index(dim-1);
                        continue;
                    }

                    //check if have reached the end of the shape
                    const Eigen::Index ii = Detail::index(index, _inv.k);
                    if(!_inv.k_shape.contains(ii))
                        break;

                    //compute phi for index
                    const auto phi = Detail::phi(_t, _x, phis, index, _inv);

                    for (size_t d = 0; d < dim; ++d) {
                        Index ni = index;
                        ni(d) += 1;
                        const Eigen::Index ii = Detail::index(ni, _inv.k);
                        phis.insert( {ii, phi(d)} );
                        //phis[ii] = phi(d);
                    }
                }

                bool done = false;
                for (Eigen::Index d = dim - 2; d >= 0; --d) {
                    index(d) += 1;
                    if (index(d) >=  _inv.k(d)) {
                        if (d == 0)
                            done = true;
                        else
                            index(d) = 0;
                    } else
                        break;
                }
                if (done) break;
            }

            return phis;
        } //compute

    }//Detail

    template<typename T, template<typename> class File>
    inline const Detail::Invariants<T> computeInvariants(const File<T>& _file) noexcept {
        using Vector = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1>;
        using Matrix = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>;

        Detail::Invariants<T> out;
        out.dimensions = _file.dimensions;
        out.k = _file.k_max;
        out.i_2_E_2 = std::complex<T>(0., 1.) / (2. * _file.epsilon * _file.epsilon);
        out.p = _file.p;
        out.q = _file.q;
        out.k_shape = _file.b_Ks;

        for(size_t t = 0; t < _file.timesteps; ++t){
            //phi0
            out.pre.push_back( std::pow(T(M_PI) * _file.epsilon * _file.epsilon, -T(_file.dimensions) / 4.) * std::pow(_file.Q[t].determinant(), -0.5) );
            out.P_Q_1.push_back( _file.P[t] * _file.Q[t].inverse() );

            out.i_E_2_p.push_back( (std::complex<T>(0., 1.) / _file.epsilon * _file.epsilon) * _file.p[t].transpose() );
            //phi
            out.Q_1.push_back( std::sqrt(2. / (_file.epsilon * _file.epsilon)) * _file.Q[t].inverse() );
            out.Q_1_Q_T.push_back( _file.Q[t].inverse() * _file.Q[t].conjugate() );
        }
        return out;
    }

    template<class T>
    std::complex<T> compute() noexcept {
        std::complex<T> res (0., 0.);
        for(Eigen::Index k = 0; k < file.value().Ks.size(); ++k){ 
            const Eigen::Index idx = Detail::index(file.value().Ks[k], file.value().k_max);
            assert(phis.contains(idx));
            const std::complex<T>& p = (*phis.find(idx)).second;
            res += file.value().c_0[w.curT](k) * p;
        } 
        res *= std::exp(std::complex<double>(0., 1.) * file.value().S(Eigen::Index(w.curT)));
    }

    template<class T>
    void findBoundingBox() noexcept {

    }

    void voxelsToFile() noexcept {

    }

    void fileToVoxels() noexcept {

    }

    void voxeling() noexcept {

    }

}//Hagedorn

#endif //TURBODORN_HPP
