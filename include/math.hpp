#ifndef MATH_HPP
#define MATH_HPP

#include "defines.hpp"

#include <unsupported/Eigen/CXX11/Tensor>

namespace Math {

    namespace Hagedorn {

        namespace Detail {

            template<class T, size_t Dim>
            [[nodiscard]] std::complex<T> phi_0 (
                    const Eigen::Matrix<T, Dim, 1>& _x, 
                    T _epsilon,
                    const Eigen::array<Eigen::Index, Dim> _k,
                    const Eigen::Matrix<std::complex<T>, Dim, 1>& _p,
                    const Eigen::Matrix<std::complex<T>, Dim, 1>& _q,
                    const Eigen::Matrix<std::complex<T>, Dim, Dim>& _Q, 
                    const Eigen::Matrix<std::complex<T>, Dim, Dim>& _P) noexcept {

                const auto xq = _x - _q;
                const auto e2 = _epsilon * _epsilon;
                const auto v1 = std::pow(M_PI * e2, - Dim/4.) * std::pow(_Q.determinant(), -0.5);
                const std::complex<T> v2 = (std::complex<T>(0., 1.) / (2. * e2) * xq.transpose() * _P * _Q.inverse() * xq);
                const std::complex<T> v3 = std::complex<T>(0., 1.) / e2 * _p.transpose() * xq;
                return v1 * std::exp(v2 + v3);
            };

            template<class T, size_t Dim>
            [[nodiscard]] Eigen::Matrix<std::complex<T>, Dim, 1> phi(
                    const Eigen::Tensor<T, Dim>& _phis,
                    const Eigen::array<Eigen::Index, Dim>& _k, 
                    Eigen::Index _dir,
                    const Eigen::Matrix<std::complex<T>, Dim, 1>& _x_q,
                    const Eigen::Matrix<std::complex<T>, Dim, Dim>& _Q_1_Qc) noexcept {
                
                using Index = Eigen::array<Eigen::Index, Dim>;
                using Vector = Eigen::Matrix<std::complex<T>, Dim, 1>;
                Vector res;

                Vector kp;
                for(size_t j = 0; j < Dim; ++j){
                    //early out
                    if(_k(j) - 1 < 0){
                        kp(j) = 0;
                        continue;
                    }

                    Index k_1 = _k;
                    k_1(j) -= 1;
                    kp(j) = std::sqrt(_k(j)) * phis(k_1);
                }

                const auto lhs = _x_q * phis(_k);
                const auto rhs = _Q_1_Qc * kp;
                const auto phi_1 =  lhs - rhs;
         
                return phi_1;
            };

            template<size_t sdim, class T, size_t Dim, size_t d>
            requires (sdim > 3)
            void loop(
                    Eigen::Tensor<T, Dim>& _phis, 
                    Eigen::array<Eigen::Index, Dim>& _i, 
                    const Eigen::Matrix<std::complex<T>, Dim, 1>& _x_q,
                    const Eigen::Matrix<std::complex<T>, Dim, Dim>& _Q_1_Qc, 
                    const Eigen::array<Eigen::Index, Dim>& _k,
                    const Eigen::Matrix<T, Dim, 1>& _ksj) noexcept {
                        
                for(size_t i = 0; i < _i(d); ++i){
                    _phis(d) = i;
                    loop<T, Dim, d - 1>(_phis, _i);
                }
            };

            template<size_t sdim, class T, size_t Dim, size_t d>
            requires (d == 1 && sdim > 3)
            void loop(
                    Eigen::Tensor<T, Dim>& _phis, 
                    Eigen::array<Eigen::Index, Dim>& _i, 
                    const Eigen::Matrix<std::complex<T>, Dim, 1>& _x_q,
                    const Eigen::Matrix<std::complex<T>, Dim, Dim>& _Q_1_Qc, 
                    const Eigen::array<Eigen::Index, Dim>& _k,
                    const Eigen::Matrix<T, Dim, 1>& _ksj) noexcept {

                using Index = Eigen::array<Eigen::Index, Dim>;
                for(size_t i = 0; i < _i(d); ++i){                
                    const auto phi_res = phi(_phis, _i, d, _x_q, _Q_1_Qc);
                    _i(d) = i;

                    for(size_t j = 0; j < Dim; ++j){
                        Index k_1 = _i;
                        k_1(j) += 1;
                        phis(k_1) = phi_res(j) * _ksj(j);
                    }
                }
            };

            //3-dims
            template<class T, size_t Dim, size_t d>
            requires (d == 3 && Dim == 3)
            void loop(
                    Eigen::Tensor<T, Dim>& _phis, 
                    Eigen::array<Eigen::Index, Dim>& _i, 
                    const Eigen::Matrix<std::complex<T>, Dim, 1>& _x_q,
                    const Eigen::Matrix<std::complex<T>, Dim, Dim>& _Q_1_Qc, 
                    const Eigen::array<Eigen::Index, Dim>& _k,
                    const Eigen::Matrix<T, Dim, 1>& _ksj) noexcept {

                using Index = Eigen::array<Eigen::Index, Dim>;
                for(size_t x = 0;  x < _k[0]; ++x) {  
                    _i[0] = x;
                    for(size_t y = 0;  x < _k[0]; ++x) {
                        _i[1] = y;
                        for(size_t z = 0;  x < _k[0]; ++x) {                                 
                            const auto phi_res = phi(_phis, _i, d, _x_q, _Q_1_Qc);
                            _i[2] = z;
                            for(size_t j = 0; j < Dim; ++j){
                                Index k_1 = _i;
                                k_1(j) += 1;
                                phis(k_1) = phi_res(j) * _ksj(j);
                            }
                        }
                    }
                }
            };

            //2-dims
            template<class T, size_t Dim, size_t d>
            requires (d == 2 && Dim == 2)
            void loop(
                    Eigen::Tensor<T, Dim>& _phis, 
                    Eigen::array<Eigen::Index, Dim>& _i, 
                    const Eigen::Matrix<std::complex<T>, Dim, 1>& _x_q,
                    const Eigen::Matrix<std::complex<T>, Dim, Dim>& _Q_1_Qc, 
                    const Eigen::array<Eigen::Index, Dim>& _k,
                    const Eigen::Matrix<T, Dim, 1>& _ksj) noexcept{

                using Index = Eigen::array<Eigen::Index, Dim>;
                for(size_t x = 0;  x < _k[0]; ++x) {  
                    _i[0] = x;
                    for(size_t y = 0;  x < _k[0]; ++x) {               
                        const auto phi_res = phi(_phis, _i, d, _x_q, _Q_1_Qc);
                        _i[1] = y;
                        for(size_t j = 0; j < Dim; ++j){
                            Index k_1 = _i;
                            k_1(j) += 1;
                            phis(k_1) = phi_res(j) * _ksj(j);
                        }
                    }
                }
            };

            //1-dims
            template<class T, size_t Dim, size_t d>
            requires (d == 1 && Dim == 1)
            void loop(
                    Eigen::Tensor<T, Dim>& _phis, 
                    Eigen::array<Eigen::Index, Dim>& _i, 
                    const Eigen::Matrix<std::complex<T>, Dim, 1>& _x_q,
                    const Eigen::Matrix<std::complex<T>, Dim, Dim>& _Q_1_Qc, 
                    const Eigen::array<Eigen::Index, Dim>& _k,
                    const Eigen::Matrix<T, Dim, 1>& _ksj) noexcept{

                using Index = Eigen::array<Eigen::Index, Dim>;
                for(size_t x = 0;  x < _k[0]; ++x) {                    
                    const auto phi_res = phi(_phis, _i, d, _x_q, _Q_1_Qc);
                    _i[0] = x;
                    for(size_t j = 0; j < Dim; ++j){
                        Index k_1 = _i;
                        k_1(j) += 1;
                        phis(k_1) = phi_res(j) * _ksj(j);
                    }
                }
            };

        } //Detail

        template<class T, size_t Dim>
        [[nodiscard]] inline Eigen::Tensor<std::complex<T>, Dim> compute(    
                const Eigen::Matrix<T, Dim, 1>& _x, 
                T _epsilon,
                const Eigen::array<Eigen::Index, Dim> _k,
                const Eigen::Matrix<std::complex<T>, Dim, 1>& _p,
                const Eigen::Matrix<std::complex<T>, Dim, 1>& _q,
                const Eigen::Matrix<std::complex<T>, Dim, Dim>& _Q, 
                const Eigen::Matrix<std::complex<T>, Dim, Dim>& _P) noexcept {

            using Index = Eigen::array<Eigen::Index, Dim>;
            using Vector = Eigen::Matrix<std::complex<T>, Dim, 1>;
            using Matrix = Eigen::Matrix<std::complex<T>, Dim, Dim>;

            Eigen::Tensor<std::complex<T>, Dim> phis(_k);

            Eigen::Matrix<T, Dim, 1> skjk;
            for(size_t i = 0; i < Dim; ++i)
                skjk(i) = sqrt(_k(i) + 1);

            const Matrix Qi = _Q.inverse();
            const Vector x_q = std::sqrt(2. / (_epsilon * _epsilon)) * Qi * (_x - _q);
            const Matrix Q_1_Qc = Qi * _Q.conjugate();

            //loop over ks
            Index i;
            i.fill(0);
            loop<T, Dim, Dim>(phis, i, x_q, Q_1_Qc, _k, skjk);

            return phis;
        };

        template<class T>
        [[nodiscard]] Eigen::Matrix<T, 3, 1> C_to_HSV(const std::complex<T>& _c){
            using Vector = Eigen::Matrix<T, 3, 1>;
            const T phase = std::arg(_c);
            Vector hsv (0.5 * std::fmod(phase + 2. * M_PI, 2. * M_PI) / M_PI, 1., 1.);
            const std::complex<T> modulus = std::abs(_c);
            //lightness
            hsv(2) = 2. * std::atan2(modulus.real, 1.) / M_PI;
            //saturation
            const T l = hsv[2];
            hsv[1] = (l <= 0.5) ? 2 * l : 2. * (1. - l);
            return hsv;
        };

    } //Hagedorn

    namespace Intersector {

        template <class T>
        struct Ray {
            Ray(const glm::vec<3, T, glm::defaultp>& _orig, const glm::vec<3, T, glm::defaultp>& _dir) : 
                orig(_orig), dir(_dir) {}
            glm::vec<3, T, glm::defaultp> orig;
            glm::vec<3, T, glm::defaultp> dir;
        };

        template <class T>
        struct Sphere {
            Sphere(T _radius, const glm::vec<3, T, glm::defaultp>& _centre, const glm::vec<3, T, glm::defaultp>& _col) :
                radius(_radius), centre(_centre), color(_col) {}
            T radius;
            glm::vec<3, T, glm::defaultp> centre;
            glm::vec<3, T, glm::defaultp> color;
        };

        template <class T>
        struct Plane {
            Plane(const glm::vec<3, T, glm::defaultp>& _p0, const glm::vec<3, T, glm::defaultp>& _n, T _hw, T _hh, const glm::vec<3, T, glm::defaultp>& _col) :
                n(_n), p0(_p0), hh(_hh), hw(_hw), color(_col) {}
            glm::vec<3, T, glm::defaultp> n;
            glm::vec<3, T, glm::defaultp> p0;  //assumed the centre of the plane
            T hh, hw;
            glm::vec<3, T, glm::defaultp> color;
        };

        template <class T>
        struct Zylinder {
            Zylinder (T _radius, const glm::vec<3, T, glm::defaultp>& _p1, const glm::vec<3, T, glm::defaultp>& _p2, const glm::vec<3, T, glm::defaultp>& _col) :
                radius(_radius), p1(_p1), p2(_p2), color(_col) {}
            T radius;
            glm::vec<3, T, glm::defaultp> p1;
            glm::vec<3, T, glm::defaultp> p2;
            glm::vec<3, T, glm::defaultp> color;
        };

        template <class T>
        struct Wheel {
            Wheel (const glm::vec<3, T, glm::defaultp>& _centre, T _r, T _R, T _thickness, const glm::mat<3, 3, T, glm::defaultp>& _rot, const glm::vec<3, T, glm::defaultp>& _col) :
                centre(_centre), r(_r), R(_R), thickness(_thickness), rot(_rot), col(_col){}
            glm::vec<3, T, glm::defaultp> centre;
            T r, R, thickness;
            glm::mat<3, 3, T, glm::defaultp> rot;
            glm::vec<3, T, glm::defaultp> col;
        };

        template <class T>
        struct Circle {
            Circle(const glm::vec<3, T, glm::defaultp>& _centre, T _radius, glm::vec<3, T, glm::defaultp> _col) :
                centre(_centre), radius(_radius), color(_col) {}
            glm::vec<3, T, glm::defaultp> centre;
            T radius;
            glm::vec<3, T, glm::defaultp> color;
        };

        template <class T>
        [[nodiscard]] bool sphere(const Sphere<T>& _s, const Ray<T>& _ray, glm::vec<3, T, glm::defaultp>& _intersect, glm::vec<3, T, glm::defaultp>& _normal, T& _t) noexcept {
            using vec3 = glm::vec<3, T, glm::defaultp>;
            const vec3 oc = _ray.orig - _s.centre;
            const T d = glm::dot(oc, _ray.dir);
            const T c = (glm::dot(oc, oc) - _s.radius * _s.radius);
            if (c > T(0.) && d > T(0)) return false;
            const T delta = d * d - c;
            if(delta < T(0.)) return false;
            _t = -d - sqrt(delta);
            if (_t < T(0.)) _t = T(0);
            _intersect = _ray.orig + _t*_ray.dir;
            _normal = glm::normalize(_intersect - _s.centre);
            return delta >= T(0.);
        };

        template <class T>
        [[nodiscard]] bool plane(const Plane<T>& _p, const Ray<T>& _ray, glm::vec<3, T, glm::defaultp>& _intersect, glm::vec<3, T, glm::defaultp>& _normal, T& _t) noexcept {
            const T s = glm::dot(_p.n, glm::vec<3, T, glm::defaultp>(_p.p0 - _ray.orig));
            _t = glm::dot(_ray.dir, _p.n);
            const T d = _t != T(0) ? s / _t : T(0.);
            _intersect = _ray.orig + d * _ray.dir;
            _normal = _p.n;
            const glm::vec<3, T, glm::defaultp> dist = _intersect - _p.p0;
            return d != T(0.) && glm::abs(dist.x) <= _p.hw && glm::abs(dist.y) <= _p.hh;
        };

        template <class T>
        [[nodiscard]] bool zylinder(const Zylinder<T>& _p, const Ray<T>& _ray, glm::vec<3, T, glm::defaultp>& _intersect, glm::vec<3, T, glm::defaultp>& _normal, T& _t) noexcept {
            using vec3 = glm::vec<3, T, glm::defaultp>;
            const vec3 ba = _p.p2 - _p.p1;
            const vec3  oc = _ray.orig - _p.p1;

            const T baba = glm::dot(ba,ba);
            const T bard = glm::dot(ba,_ray.dir);
            const T baoc = glm::dot(ba,oc);
            
            const T k2 = baba - bard * bard;
            const T k1 = baba * glm::dot(oc,_ray.dir) - baoc * bard;
            const T k0 = baba * glm::dot(oc,oc) - baoc * baoc - _p.radius * _p.radius * baba;
            const T h = k1 * k1 - k2 * k0;

            if(h < T(0.)) return false; //e1
            h = glm::sqrt(h);
            _t = (-k1 - h) / k2;

            if(_t < T(0)) return false;//e2

            // body
            const T y = baoc + _t * bard;
            if(y > T(0.) && y < baba){ //e3
                _intersect = _ray.orig + _t * _ray.dir;
                _normal = glm::normalize((oc + _t * _ray.dir - ba * y / baba) / _p.radius);
                return true;
            }

            // caps
            _t = ( ( (y < T(0.)) ? T(0.) : baba) - baoc)/bard;
            if( glm::abs(k1 + k2 * _t) < h){ //e4
                _intersect = _ray.orig + _t * _ray.dir;
                _normal = glm::normalize(ba * glm::sign(y) / baba);
                return true;
            }

            return false; //no intersection
        };

        template <class T>
        [[nodiscard]] bool circle(const Circle<T>& _p, const Ray<T>& _ray, glm::vec<3, T, glm::defaultp>& _intersect, glm::vec<3, T, glm::defaultp>& _normal, T& _t) noexcept {
            using vec3 = glm::vec<3, T, glm::defaultp>;
            const T s = glm::dot(_p.n, vec3(_p.p0 - _ray.orig));
            _t = glm::dot(_ray.dir, _p.n);
            const T d = _t != T(0.) ? s / _t : T(0.);
            _intersect = _ray.orig + d * _ray.dir;
            _normal = _p.n;
            const vec3 dist = _intersect - _p.p0;
            const T dist2 = glm::dot(dist, dist);
            return d != T(0.) && dist2 >= _p.radius_inner*_p.radius_inner && dist2 <= _p.radius_outer*_p.radius_outer;
        };

        namespace SDF {

            template <class T>
            [[nodiscard]] T wheel(const Wheel<T>& _wheel, const glm::vec<3, T, glm::defaultp>& _pos) noexcept {
                using vec2 = glm::vec<2, T, glm::defaultp>;
                using vec3 = glm::vec<3, T, glm::defaultp>;
                using vec4 = glm::vec<4, T, glm::defaultp>;
                const vec3 pos = _wheel.rot * _pos;
                const vec3 q = glm::max(glm::abs(pos) - vec3(T(0.), _wheel.thickness, T(0.)), T(0.));
                const vec3 v = vec4( glm::max(q, T(0.)), glm::min( glm::max(q.x, glm::max(q.y, q.z)), T(0.) ) );
                const vec3 C = vec3(_wheel.cX, _wheel.cY, _wheel.cZ);
                const vec2 t = vec2(glm::length(v.xz - C.xz) - _wheel.R, v.y - C.y);
                return glm::length(t) - _wheel.r;
            };

            template <class T>
            [[nodiscard]] T sphere(const Sphere<T>& _sphere, const glm::vec<3, T, glm::defaultp>& _pos) noexcept {
                using vec3 = glm::vec<3, T, glm::defaultp>;
                const vec3 sc = vec3(_sphere.centreX, _sphere.centreY, _sphere.centreZ);
                return glm::length(_pos - sc) - _sphere.radius;
            };

            template <class T>
            [[nodiscard]] T line(const Zylinder<T>& _p, const glm::vec<3, T, glm::defaultp>& _pos) {
                using vec3 = glm::vec<3, T, glm::defaultp>;
                const vec3 p1 = vec3(_p.p1X, _p.p1Y, _p.p1Z);
                const vec3 p2 = vec3(_p.p2X, _p.p2Y, _p.p2Z);
                const vec3  ba = p2 - p1;
                const vec3  pa = _pos - p1;
                const T h = glm::clamp( glm::dot(pa,ba)/glm::dot(ba,ba), T(0.), T(1.) );
                return glm::length( pa - ba*h ) - _p.radius;
            };

            template <class T>
            [[nodiscard]] T zylinder(const Zylinder<T>& _p, const glm::vec<3, T, glm::defaultp>& _pos) {
                using vec3 = glm::vec<3, T, glm::defaultp>;
                const vec3 p1 = vec3(_p.p1X, _p.p1Y, _p.p1Z);
                const vec3 p2 = vec3(_p.p2X, _p.p2Y, _p.p2Z);
                const vec3 ba = p2 - p1;
                const vec3 pa = _pos - p1;
                const T baba = glm::dot(ba, ba);
                const T paba = glm::dot(pa, ba);
                const T x = glm::length(pa * baba - ba * paba) - _p.radius * baba;
                const T y = glm::abs(paba-baba * T(0.5)) - baba * T(0.5);
                const T x2 = x * x;
                const T y2 = y * y * baba;   
                const T d = (glm::max(x, y) < T(0.)) ? -glm::min(x2, y2) : ( ((x > T(0.)) ? x2 : T(0.)) + ((y > T(0.)) ? y2 : T(0.)) );   
                return glm::sign(d) * glm::sqrt(glm::abs(d)) / baba;
            };

        }
    }

}

#endif /* MATH_HPP */
