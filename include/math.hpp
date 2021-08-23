#ifndef MATH_HPP
#define MATH_HPP

#include "defines.hpp"

namespace Math {

    namespace Hagedorn {

        namespace Detail {

            template<class T>
            struct Invariants {
                size_t dimensions;
                Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1> k; //k extends
                std::unordered_map<Eigen::Index, bool> k_shape; //lookup for shape
                std::vector<Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1>> p;
                std::vector<Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1>> q;
                //phi_0
                std::vector<std::complex<T>> pre;
                std::complex<T> i_2_E_2;
                std::vector<Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>> P_Q_1;
                std::vector<Eigen::Matrix<std::complex<T>, 1, Eigen::Dynamic>> i_E_2_p;
                //phi
                //sqrt*Q-1
                std::vector<Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>> Q_1;
                //Q-1*QT
                std::vector<Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>> Q_1_Q_T;
            };

            // [x_n, ..., x_2, x_1] -> N
            [[nodiscard]] Eigen::Index index(const Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1>& /*_i*/, const Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1>& /*_e*/) noexcept;

            template<class T>
            [[nodiscard]] std::complex<T> phi_0 ( size_t /*_t*/, const Eigen::Matrix<T, Eigen::Dynamic, 1>& /*_x*/, const Detail::Invariants<T>& /*_inv*/) noexcept;

            template<class T>
            [[nodiscard]] Eigen::Matrix<std::complex<T>, -1, 1> phi (
                    size_t /*_t*/,
                    const Eigen::Matrix<T, Eigen::Dynamic, 1>& /*_x*/,
                    const std::unordered_map<Eigen::Index, std::complex<T>>& /*_phis*/,
                    const Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1>& /*_index*/,
                    const Detail::Invariants<T>& /*_inv*/
            ) noexcept; 

        } //Detail

        template<typename T, template<typename> class File>
        [[nodiscard]] const Detail::Invariants<T> computeInvariants(const File<T>& /*_file*/) noexcept;

        template<class T>
        [[nodiscard]] std::unordered_map<Eigen::Index, std::complex<T>> compute ( size_t /*_t*/, const Eigen::Matrix<T, Eigen::Dynamic, 1>& /*_x*/, const Detail::Invariants<T>& /*_inv*/) noexcept;

    } //Hagedorn

    template<class Vec, class T>
    [[nodiscard]] bool intersect(const Vec& /*_r_o*/, const Vec& /*_r_o*/, const Vec& /*_low*/, const Vec& /*_high*/, T /*_tMax*/, T& /*_t*/) noexcept;

}

template<typename T, template<typename> class File>
inline const Math::Hagedorn::Detail::Invariants<T> Math::Hagedorn::computeInvariants(const File<T>& _file) noexcept {
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

template <class T>
inline std::unordered_map<Eigen::Index, std::complex<T>> Math::Hagedorn::compute (size_t _t, const Eigen::Matrix<T, Eigen::Dynamic, 1>& _x, const Detail::Invariants<T>& _inv) noexcept {

    using Index = Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1>;
    using Vector = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1>;
    using Matrix = Eigen::Matrix<std::complex<T>, Eigen::Dynamic, Eigen::Dynamic>;

    const size_t dim = _inv.dimensions;

    size_t size = _inv.k(0) + 1;
    for (size_t i = 1; i < _inv.k.size(); ++i)
        size *= (_inv.k(i)+1);

    std::unordered_map<Eigen::Index, std::complex<T>> phis;
    //phis.resize(size);

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
}; //compute

/*
    i: index
    e: extends, # units
*/
inline Eigen::Index Math::Hagedorn::Detail::index(const Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1>& _i, const Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1>& _e) noexcept {
    assert(_i.size() == _e.size());
    Eigen::Index out = _i(0);
    for (Eigen::Index k = 1; k < _i.size(); ++k) {
        out *= _e(k);
        out += _i(k);
    }
    return out;
};  //index

template <class T>
inline std::complex<T> Math::Hagedorn::Detail::phi_0 (size_t _t, const Eigen::Matrix<T, Eigen::Dynamic, 1>& _x, const Detail::Invariants<T>& _inv) noexcept {
    const Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> xq = _x - _inv.q[_t];
    const Eigen::Matrix<std::complex<T>, 1, Eigen::Dynamic> xqt = xq.transpose();
    const std::complex<T> e1 = _inv.i_2_E_2 * xqt * _inv.P_Q_1[_t] * xq;
    const std::complex<T> e2 = _inv.i_E_2_p[_t] * xq;
    return _inv.pre[_t] * std::exp(e1 + e2);
};  //phi_0

template <class T>
inline Eigen::Matrix<std::complex<T>, Eigen::Dynamic, 1> Math::Hagedorn::Detail::phi (
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
        //early out
        if (_index(j) - 1 < 0) {
            kp(j) = { 0., 0. };
            continue;
        }

        Index k_1 = _index;
        k_1(j) = std::max(k_1(j) - 1, 0ll);
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
};  //phi

template<class Vec, class T>
inline bool Math::intersect(const Vec& _r_o, const Vec& _r_d, const Vec& _low, const Vec& _high, T _tmax, T& _t) noexcept {
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
}; //intersect

#endif /* MATH_HPP */
