#ifndef MATH_HPP
#define MATH_HPP

#include <Eigen/Eigen>
#include <Eigen/Core>

#include <vector>
#include <memory>
#include <functional>
#include <iostream>
#include <complex>
#include <algorithm>

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

                Invariants() noexcept = default;
                Invariants(Invariants&& _in) noexcept = default;               
                [[nodiscard]] Invariants& operator=(Invariants&& _in) noexcept = default;
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

    template <class T>
    [[nodiscard]] Eigen::Matrix<T, 3, 1> c_to_HSL(T /*_max*/, const std::complex<T>& /*_c*/) noexcept;

    template <class T>
    [[nodiscard]] Eigen::Matrix<T, 3, 1> HSL_to_RGB_rad(const Eigen::Matrix<T, 3, 1>& /*_hsl*/) noexcept;

    template <class T>
    [[nodiscard]] Eigen::Matrix<T, 3, 1> HSL_to_RGB_deg(const Eigen::Matrix<T, 3, 1>& /*_hsl*/) noexcept;

    template <class T>
    [[nodiscard]] Eigen::Matrix<T, 3, 1> rgb_to_gs(const Eigen::Matrix<T, 3, 1>& /*_rgb*/) noexcept;

    template <class T>
    [[nodiscard]] constexpr T depth(const Eigen::Array<T, 4, 1>& /*_ps*/) noexcept;

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
inline std::unordered_map<Eigen::Index, std::complex<T>> Math::Hagedorn::compute (
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
};//Math::Hagedorn::compute

/*
    i: index
    e: extends, # units
*/
inline Eigen::Index Math::Hagedorn::Detail::index(
    const Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1>& _i, 
    const Eigen::Matrix<Eigen::Index, Eigen::Dynamic, 1>& _e
) noexcept {
    assert(_i.size() == _e.size());
    Eigen::Index out = _i(0);
    for (Eigen::Index k = 1; k < _i.size(); ++k) {
        out *= _e(k);
        out += _i(k);
    }
    return out;
}; //Math::Hagedorn::Detail::index







template <class T>
inline constexpr T Math::depth(const Eigen::Array<T, 4, 1>& _p) noexcept{
    return (_p.z / _p.z + 1.) * 0.5;
} //depth

#endif /* MATH_HPP */
