#ifndef MATH_HPP
#define MATH_HPP

#include "defines.hpp"

namespace Math {

    namespace Hagedorn {

        namespace Detail {

            // [x_n, ..., x_2, x_1] -> N
            [[nodiscard]] Eigen::Index index(const Eigen::Matrix<Eigen::Index, -1, 1>& _i, const Eigen::Matrix<Eigen::Index, -1, 1>& _e) noexcept;

            template<class T>
            [[nodiscard]] std::complex<T> phi_0 (
                    const Eigen::Matrix<T, -1, 1>& _x, 
                    T _epsilon,
                    const Eigen::Matrix<std::complex<T>, -1, 1>& _p,
                    const Eigen::Matrix<std::complex<T>, -1, 1>& _q,
                    const Eigen::Matrix<std::complex<T>, -1, -1>& _Q, 
                    const Eigen::Matrix<std::complex<T>, -1, -1>& _P) noexcept;

            template<class T>
            [[nodiscard]] Eigen::Matrix<std::complex<T>, -1, 1> phi (
                    const std::vector<std::complex<T>>& _phis,
                    const Eigen::Matrix<Eigen::Index, -1, 1>& _k, 
                    const Eigen::Matrix<Eigen::Index, -1, 1>& _index,
                    const Eigen::Matrix<std::complex<T>, -1, 1>& _x_q,
                    const Eigen::Matrix<std::complex<T>, -1, -1>& _Q_1_Qc) noexcept; 

        } //Detail

        template<class T>
        [[nodiscard]] std::vector<std::complex<T>> compute (
            const Eigen::Matrix<T, -1, 1>& _x, 
            T _epsilon,
            const Eigen::Matrix<Eigen::Index, -1, 1>& _k,
            const Eigen::Matrix<std::complex<T>, -1, 1>& _p,
            const Eigen::Matrix<std::complex<T>, -1, 1>& _q,
            const Eigen::Matrix<std::complex<T>, -1, -1>& _Q, 
            const Eigen::Matrix<std::complex<T>, -1, -1>& _P) noexcept;


    } //Hagedorn

    template<class Vec, class T>
    [[nodiscard]] bool intersect(const Vec& /*_r_o*/, const Vec& /*_r_o*/, const Vec& /*_low*/, const Vec& /*_high*/, T /*_tMax*/, T& /*_t*/) noexcept;

}

template <class T>
inline std::vector<std::complex<T>> Math::Hagedorn::compute (
    const Eigen::Matrix<T, -1, 1>& _x,
    T _epsilon,
    const Eigen::Matrix<Eigen::Index, -1, 1>& _k,
    const Eigen::Matrix<std::complex<T>, -1, 1>& _p,
    const Eigen::Matrix<std::complex<T>, -1, 1>& _q,
    const Eigen::Matrix<std::complex<T>, -1, -1>& _Q,
    const Eigen::Matrix<std::complex<T>, -1, -1>& _P) noexcept {

    using Index = Eigen::Matrix<Eigen::Index, -1, 1>;
    using Vector = Eigen::Matrix<std::complex<T>, -1, 1>;
    using Matrix = Eigen::Matrix<std::complex<T>, -1, -1>;

    const size_t dim = _k.rows();

    size_t size = _k(0) + 1;
    for (size_t i = 1; i < _k.size(); ++i)
        size *= (_k(i)+1);

    std::vector<std::complex<T>> phis;
    phis.resize(size);

    const Matrix Qi = _Q.inverse();
    const Vector x_q = std::sqrt(2. / (_epsilon * _epsilon)) * Qi * (_x - _q);
    const Matrix Q_1_Qc = Qi * _Q.conjugate();

    //iterate over ks
    Index index(dim);
    index.fill(0);

    bool first = true;

    while (true) {
        for (index(dim-1) = 0; index(dim-1) <= _k(dim-1); ++index(dim-1)) { 
            if (first) {
                first = false;
                const auto phi0 = Detail::phi_0(_x, _epsilon, _p, _q, _Q, _P);
                phis[0] = phi0;
                continue;
            }

            const auto phi = Detail::phi(phis, _k, index, x_q, Q_1_Qc);
            //std::cout << phi << std::endl;

            for (size_t d = 0; d < dim; ++d) {
                Index ni = index;
                ni(d) += 1;
                const Eigen::Index ii = Detail::index(ni, _k);
                phis[ii] = phi(d);
            }
        }

        bool done = false;
        for (Eigen::Index d = dim - 2; d >= 0; --d) {
            index(d) += 1;
            if (index(d) >= _k(d)) {
                if (d == 0)
                    done = true;
                else
                    index(d) = 0;
            } else
                break;
        }
        if (done) break;
    }

    //for (size_t i = 0; i < dim; ++i) {
        //const std::complex<T> ki(T(_k(i)) + 1., 0.);
      //  phis(i) /= sqrt(ki);
    //}

    return phis;
}; //compute

/*
    i: index
    e: extends, # units
*/
inline Eigen::Index Math::Hagedorn::Detail::index(const Eigen::Matrix<Eigen::Index, -1, 1>& _i, const Eigen::Matrix<Eigen::Index, -1, 1>& _e) noexcept {
    assert(_i.size() == _e.size());
    Eigen::Index out = _i(0);
    for (Eigen::Index k = 1; k < _i.size(); ++k) {
        out *= _e(k);
        out += _i(k);
    }
    return out;
};  //index

template <class T>
inline std::complex<T> Math::Hagedorn::Detail::phi_0 (
    const Eigen::Matrix<T, -1, 1>& _x,
    T _epsilon,
    const Eigen::Matrix<std::complex<T>, -1, 1>& _p,
    const Eigen::Matrix<std::complex<T>, -1, 1>& _q,
    const Eigen::Matrix<std::complex<T>, -1, -1>& _Q,
    const Eigen::Matrix<std::complex<T>, -1, -1>& _P) noexcept {

    const size_t dims = _x.rows();
    const auto xq = _x - _q;
    const T e2 = _epsilon * _epsilon;

    const auto v1_a = std::pow(_Q.determinant(), -0.5);
    const auto v1_b = std::pow(T(M_PI) * e2, -T(dims) / 4.);

    const std::complex<T> v2 = ((std::complex<T>(0., 1.) / (2. * e2)) * xq.transpose() * _P * _Q.inverse() * xq);
    const std::complex<T> v3 = (std::complex<T>(0., 1.) / e2) * _p.transpose() * xq;

    return v1_a * v1_b * std::exp(v2 + v3);
};  //phi_0

template <class T>
inline Eigen::Matrix<std::complex<T>, -1, 1> Math::Hagedorn::Detail::phi (
    const std::vector<std::complex<T>>& _phis,
    const Eigen::Matrix<Eigen::Index, -1, 1>& _k,
    const Eigen::Matrix<Eigen::Index, -1, 1>& _index,
    const Eigen::Matrix<std::complex<T>, -1, 1>& _x_q,
    const Eigen::Matrix<std::complex<T>, -1, -1>& _Q_1_Qc) noexcept {

    using Index = Eigen::Matrix<Eigen::Index, -1, 1>;
    using Vector = Eigen::Matrix<std::complex<T>, -1, 1>;
    using Matrix = Eigen::Matrix<std::complex<T>, -1, -1>;

    const size_t dim = _k.rows();

    //std::cout << _index << std::endl << std::endl << _k << std::endl;

    Vector res (dim);

    Vector kp (dim);
    for (size_t j = 0; j < dim; ++j) {
        //early out
        if (_index(j) - 1 < 0) {
            kp(j) = 0;
            continue;
        }

        Index k_1 = _index;
        k_1(j) = std::max(k_1(j) - 1, 0ll);
        //std::cout << _index << std::endl << k_1 << std::endl;
        //std::cout << _k << std::endl << std::endl;
        const Eigen::Index ii = Detail::index(k_1, _k);
        kp(j) = std::sqrt(_index(j)) * _phis[ii];  
    }

    const Eigen::Index ii = Detail::index(_index, _k);
    const auto lhs = _x_q *_phis[ii];
    const auto rhs = _Q_1_Qc * kp;
    auto phi_t = lhs - rhs;

    Vector phi (dim);
    for(size_t i = 0; i < dim; ++i){
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
