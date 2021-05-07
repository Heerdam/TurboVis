#ifndef MATH_HPP
#define MATH_HPP

#include "defines.hpp"

namespace Math {

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
            T line(const Zylinder<T>& _p, const glm::vec<3, T, glm::defaultp>& _pos) {
                using vec3 = glm::vec<3, T, glm::defaultp>;
                const vec3 p1 = vec3(_p.p1X, _p.p1Y, _p.p1Z);
                const vec3 p2 = vec3(_p.p2X, _p.p2Y, _p.p2Z);
                const vec3  ba = p2 - p1;
                const vec3  pa = _pos - p1;
                const T h = glm::clamp( glm::dot(pa,ba)/glm::dot(ba,ba), T(0.), T(1.) );
                return glm::length( pa - ba*h ) - _p.radius;
            };

            template <class T>
            T zylinder(const Zylinder<T>& _p, const glm::vec<3, T, glm::defaultp>& _pos) {
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
