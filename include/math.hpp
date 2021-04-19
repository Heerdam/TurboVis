#ifndef MATH_HPP
#define MATH_HPP

#include "defines.hpp"

namespace Math {

    namespace Intersector {

        template <class Vector>
        struct Ray {
            Ray(const Vector& _orig, const Vector& _dir) : 
                orig(_orig), dir(_dir) {}
            Vector orig;
            Vector dir;
        };

        template <class Vector, class T>
        struct Sphere {
            Sphere(T _radius, const Vector& _centre) :
                radius(_radius), centre(_centre) {}
            T radius;
            Vector centre;
        };

        template <class Vector, class T>
        struct Plane {
            Plane(const Vector& _p0, const Vector& _n, T _hw, T _hh) :
                n(_n), p0(_p0), hh(_hh), hw(_hw) {}
            Vector n;
            Vector p0;  //assumed the centre of the plane
            T hh, hw;
        };

        template <class Vector, class T>
        struct Zylinder {
            Zylinder (T _radius, const Vector& _p1, const Vector& _p2) :
                radius(_radius), p1(_p1), p2(_p2) {}
            T radius;
            Vector p1;
            Vector p2;
        };

        template <class Vector, class T>
        struct Circle {
            Circle(const Vector& _p0, const Vector& radius_inner, T _hw, T radius_outer) :
                n(_n), p0(_p0), radius_inner(_radius_inner), radius_outer(_radius_outer) {}
            Vector n;
            Vector p0;  //assumed the centre of the plane
            T radius_outer, radius_inner;
        };

        template <class Vector, class T>
        [[nodiscard]] bool sphere(const Sphere<Vector, T>& _s, const Ray<Vector>& _ray, Vector& _intersect, Vector& _normal, T& _t) noexcept {
            const Vector oc = _ray.orig - _s.centre;
            const T d = glm::dot(oc, _ray.dir);
            const T c = (glm::dot(oc, oc) - _s.radius * _s.radius);
            if (c > T(0) && d > T(0)) return false;
            const T delta = d * d - c;
            if(delta < T(0)) return false;
            _t = -d - sqrt(delta);
            if (_t < T(0)) _t = T(0);
            _intersect = _ray.orig + _t*_ray.dir;
            _normal = glm::normalize(Vector(_intersect - _s.centre));
            return delta >= T(0);
        };

        template <class Vector, class T>
        [[nodiscard]] bool plane(const Plane<Vector, T>& _p, const Ray<Vector>& _ray, Vector& _intersect, Vector& _normal, T& _t) noexcept {
            const T s = glm::dot(_p.n, Vector(_p.p0 - _ray.orig));
            _t = glm::dot(_ray.dir, _p.n);
            const T d = _t != T(0) ? s / _t : T(0);
            _intersect = _ray.orig + d * _ray.dir;
            _normal = _p.n;
            const Vector dist = _intersect - _p.p0;
            return d != T(0) && glm::abs(dist.x) <= _p.hw && glm::abs(dist.y) <= _p.hh;
        };

        template <class Vector, class T>
        [[nodiscard]] bool zylinder(const Zylinder<Vector, T>& _p, const Ray<Vector>& _ray, Vector& _intersect, Vector& _normal, T& _t) noexcept {
            const Vector ba = _p.p2 - _p.p1;
            const Vector  oc = _ray.orig - _p.p1;

            const T baba = glm::dot(ba,ba);
            const T bard = glm::dot(ba,_ray.dir);
            const T baoc = glm::dot(ba,oc);
            
            const T k2 = baba - bard * bard;
            const T k1 = baba * glm::dot(oc,_ray.dir) - baoc * bard;
            const T k0 = baba * glm::dot(oc,oc) - baoc * baoc - _p.radius * _p.radius * baba;
            const T h = k1 * k1 - k2 * k0;

            if(h < T(0)) return false; //e1
            h = glm::sqrt(h);
            _t = (-k1 - h) / k2;

            if(_t < T(0)) return false;//e2

            // body
            const T y = baoc + _t * bard;
            if(y > T(0) && y < baba){ //e3
                _intersect = _ray.orig + _t * _ray.dir;
                _normal = glm::normalize((oc + _t * _ray.dir - ba * y / baba) / _p.radius);
                return true;
            }

            // caps
            _t = ( ( (y < T(0)) ? T(0) : baba) - baoc)/bard;
            if( glm::abs(k1 + k2 * _t) < h){ //e4
                _intersect = _ray.orig + _t * _ray.dir;
                _normal = glm::normalize(ba * glm::sign(y) / baba);
                return true;
            }

            return false; //no intersection
        };

        template <class Vector, class T>
        [[nodiscard]] bool circle(const Circle<Vector, T>& _p, const Ray<Vector>& _ray, Vector& _intersect, Vector& _normal, T& _t) noexcept {
            const T s = glm::dot(_p.n, Vector(_p.p0 - _ray.orig));
            _t = glm::dot(_ray.dir, _p.n);
            const T d = _t != T(0) ? s / _t : T(0);
            _intersect = _ray.orig + d * _ray.dir;
            _normal = _p.n;
            const Vector dist = _intersect - _p.p0;
            const T dist2 = glm::dot(dist, dist);
            return d != T(0) && dist2 >= _p.radius_inner*_p.radius_inner && dist2 <= _p.radius_outer*_p.radius_outer;
        };

        namespace Branchless {

        }
    }

}

#endif /* MATH_HPP */
