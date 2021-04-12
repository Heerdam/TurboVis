
#include "../../include/defines.hpp"

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>

#include <lodepng.h>

using vec3 = glm::vec3;
using vec2 = glm::vec2;

using mat4 = glm::mat4;

using namespace glm;

#define EPSILON 0.00001

template<class T>
[[nodiscard]] constexpr T sqr(T _v) noexcept{
    return _v * _v;
};

template<class Vec>
[[nodiscard]] constexpr float dist2(const Vec& _v) noexcept{
    return dot(_v, _v);
};

[[nodiscard]] vec3 operator*(const mat4& _cam, const vec3& _vec) noexcept {
    const vec4 t = _cam * vec4(_vec.x, _vec.y, _vec.z, 1.f);
    return vec3(t.x, t.y, t.z);
};

struct Material {
    Material(const vec3& _c) : col(_c) {}
    vec3 col;
};

struct DirLight {
    vec3 dir;
    vec3 col;
    float strength;
};

struct Ambiente {
    vec3 col;
    float strength;
};

class Camera {

    vec3 dir;
    vec3 up;
    vec3 right;

    const vec3 upAxis = vec3(0.f, 1.f, 0.f);

public:

    vec3 pos;
    vec3 target;

    const float near;
    const float far;
    const float fovX;
    const float fovY;
    const int64_t w;
    const int64_t h;

    void update() noexcept {
        dir = normalize(target - pos);
        right = normalize(cross(dir, upAxis));
        up = normalize(cross(right, dir));
    };

    [[nodiscard]] Vec3& operator*(Vec3& _in) const noexcept {
        _in = normalize(_in.x*right + _in.y*up + _in.z*dir);
        return _in;
    };

    Camera(size_t _w, size_t _h, float _fov, float _near, float _far) noexcept
        : w(_w), h(_h), fovX(_fov), fovY(_w / _h * _fov), near(_near), far(_far) {}

};

struct Ray {
    vec3 orig;
    vec3 dir;
    vec3 end;
};

struct Sphere {
    Sphere(float _r, const vec3& _c, const Material& _m) : 
        radius(_r), centre(_c), mat(_m) {}

    float radius;
    vec3 centre;
    Material mat;
};

struct Plane {
    Plane(const vec3& _n, const vec3& _p0, float _hh, float _hw, const Material& _mat) :
        n(_n), p0(_p0), hh(_hh), hw(_hw), mat(_mat) {}

    vec3 n;
    vec3 p0; //assumed the centre of the plane
    float hh, hw;
    Material mat;
};

struct Zylinder {
    Zylinder(float _r, const vec3& _p1, const vec3& _p2, const Material& _mat) :
        radius(_r), p1(_p1), p2(_p2), mat(_mat) {}

    float radius;
    vec3 p1;
    vec3 p2;
    Material mat;
};

bool sphere(const Sphere& _s, const Ray& _ray, vec3& _intersect, vec3& _normal, float& _t) {
    const vec3 oc = _ray.orig - _s.centre;
    const float d = dot(oc, _ray.dir);
    const float c = (dot(oc, oc) - _s.radius * _s.radius);
    if (c > 0.f && d > 0.f) return false;
    const float delta = d * d - c;
    if(delta < 0.f) return false;
    _t = -d - sqrt(delta);
    if (_t < 0.f) _t = 0.f;
    _intersect = _ray.orig + _t*_ray.dir;
    _normal = normalize(vec3(_intersect - _s.centre));
    return delta >= 0.f;
};

bool sphere_bl(const Sphere& _s, const Ray& _ray, vec3& _intersect, vec3& _normal, float& _t) {
    const vec3 oc = _ray.orig - _s.centre;
    const float d = dot(oc, _ray.dir);
    const float c = (dot(oc, oc) - _s.radius * _s.radius);
    const int e1 = int(c > 0.f && d > 0.f);
    const float delta = d * d - c;
    const int e2 = int(delta < 0.f);
    _t = -d - sqrt(delta);
    const int e3 = int((_t < 0.f));
    _t = (1 - e3) * _t;
    _intersect = _ray.orig + _t*_ray.dir;
    _normal = normalize(vec3(_intersect - _s.centre));
    const int e4 = int(delta >= 0.f);
    const int out = (1 - e1)*(1 - e2)*e4;
    return bool(out);
};

bool plane(const Plane& _p, const Ray& _ray, vec3& _intersect, vec3& _normal, float& _t){
    const float s = dot(_p.n, vec3(_p.p0 - _ray.orig));
    _t = dot(_ray.dir, _p.n);
    const float d = _t != 0.f ? s / _t : 0.f;
    _intersect = _ray.orig + d * _ray.dir;
    _normal = _p.n;
    const vec3 dist = _intersect - _p.p0;
    return d != 0.f && abs(dist.x) <= _p.hw && abs(dist.y) <= _p.hh;
};

bool plane_bl(const Plane& _p, const Ray& _ray, vec3& _intersect, vec3& _normal, float& _t){
    const float s = dot(_p.n, vec3(_p.p0 - _ray.orig));
    _t = dot(_ray.dir, _p.n);
    const int e1 = int(_t != 0.f);
    const float d = e1 * (s / _t);
    _intersect = _ray.orig + d * _ray.dir;
    _normal = _p.n;
    const vec3 dist = _intersect - _p.p0;
    return d != 0.f && abs(dist.x) <= _p.hw && abs(dist.y) <= _p.hh;
};

bool zylinder(const Zylinder& _p, const Ray& _ray, vec3& _intersect, vec3& _normal, float& _t){   

    const vec3 ba = _p.p2 - _p.p1;
    const vec3  oc = _ray.orig - _p.p1;

    const float baba = dot(ba,ba);
    const float bard = dot(ba,_ray.dir);
    const float baoc = dot(ba,oc);
    
    const float k2 = baba - bard * bard;
    const float k1 = baba * dot(oc,_ray.dir) - baoc * bard;
    const float k0 = baba * dot(oc,oc) - baoc * baoc - _p.radius * _p.radius * baba;
    float h = k1 * k1 - k2 * k0;

    if(h < 0.f) return false; //e1
    h = sqrt(h);
    _t = (-k1 - h) / k2;

    if(_t < 0.f) return false;//e2

    // body
    const float y = baoc + _t * bard;
    if(y > 0.f && y < baba){ //e3
        _intersect = _ray.orig + _t * _ray.dir;
        _normal = normalize((oc + _t * _ray.dir - ba * y / baba) / _p.radius);
        return true;
    }

    // caps
    _t = ( ( (y < 0.f) ? 0.f : baba) - baoc)/bard;
    if( abs(k1 + k2 * _t) < h){ //e4
        _intersect = _ray.orig + _t * _ray.dir;
        _normal = normalize(ba * sign(y) / baba);
        return true;
    }

    return false; //no intersection
};

bool zylinder_bl(const Zylinder& _p, const Ray& _ray, vec3& _intersect, vec3& _normal, float& _t){   

    const vec3 ba = _p.p2 - _p.p1;
    const vec3  oc = _ray.orig - _p.p1;

    const float baba = dot(ba,ba);
    const float bard = dot(ba,_ray.dir);
    const float baoc = dot(ba,oc);
    
    const float k2 = baba - bard * bard;
    const float k1 = baba * dot(oc,_ray.dir) - baoc * bard;
    const float k0 = baba * dot(oc,oc) - baoc * baoc - _p.radius * _p.radius * baba;
    float h = k1 * k1 - k2 * k0;

    const int e1 = int(h < 0.f);
    h = sqrt(h);
    _t = (-k1 - h) / k2;
    const int e2 = int(_t < 0.f);

    // body
    const float y = baoc + _t * bard;
    const int e3 = int(y > 0.f && y < baba);

    // caps
    const int e5 = int((y < 0.f));
    _t = ((1-e5) * baba - baoc)/bard; 
    const int e4 = int(abs(k1 + k2 * _t) < h);

    _intersect = _ray.orig + _t * _ray.dir;
    _normal = e3 * (normalize((oc + _t * _ray.dir - ba * y / baba) / _p.radius)) + (1-e3) * e4 * (normalize(ba * sign(y) / baba));

    const int out = bool((1 - e1) * (1 - e2) * (e3 + e4));
    return bool(out); 
};

[[nodiscard]] Ray cameraRay(const Camera& _cam, int64_t _u, int64_t _v) noexcept {
    const float x = (static_cast<float>(2*_u - _cam.w)/static_cast<float>(_cam.w))*std::tan(_cam.fovX);
    const float y = (static_cast<float>(2*_v - _cam.h)/static_cast<float>(_cam.h))*std::tan(_cam.fovY);
    Ray ray;
    ray.dir = _cam * vec3(x, y, 1.f);
    ray.orig = _cam.pos;
    return ray;
};

template<class Col, class Depth, class Shapes, class Func, class Camera>
void render(Col& _cb, Depth& _d, const Shapes& _shapes, Func _func, const Ray& _ray, const Ambiente& _amb, const DirLight& _light, const Camera& _cam, size_t _idx, size_t _idx_d) noexcept {
    for (size_t k = 0; k < _shapes.size(); ++k) {
        vec3 inter;
        vec3 normal;
        float t = 0.f;
        if (_func(_shapes[k], _ray, inter, normal, t)) {
            const float z = t;
            if (z < _d[_idx_d]) {
                _d[_idx_d] = z;

                const float diff = std::max(dot(normal, _light.dir), 0.f);
                const vec3 diffl = diff * _light.col * _light.strength;
                const vec3 res = (_amb.col * _amb.strength + diffl) * _shapes[k].mat.col;
                const vec3 resg = glm::pow(res, vec3(1.f / 2.2f));
                _cb[_idx] = (unsigned char)(std::clamp(resg.x * 255.f, 0.f, 255.f));
                _cb[_idx + 1] = (unsigned char)(std::clamp(resg.y * 255.f, 0.f, 255.f));
                _cb[_idx + 2] = (unsigned char)(std::clamp(resg.z * 255.f, 0.f, 255.f));
                _cb[_idx + 3] = 255;
            }
        }
    }
};

int main() {

    // ------------------- CAMERA -------------------
    Camera cam (1000, 1000, glm::radians(45.f), 0.001f, 1000.f);
    cam.pos = vec3(0.f, 0.f, -500.f);
    cam.target = vec3(0.f);
    cam.update();

    // ------------------- LIGHT -------------------
    DirLight light;
    light.dir = -normalize(vec3(1.f));
    light.col = vec3(1.f);
    light.strength = 0.75f;

    Ambiente ambiente;
    ambiente.col = vec3(1.f);
    ambiente.strength = 0.4f;
    
    // ------------------- SHAPES -------------------
    std::vector<Zylinder> zyls;
    //zyls.push_back( {100.f, vec3(0.f, -200.f, 300.f), vec3(0.f, 200.f, -300.f), Material(vec3(0.9f, 0.25f, 1.f)) });

    std::vector<Sphere> spheres;
    //spheres.push_back( {250.f, vec3(0.f), Material(vec3(0.3f, 0.25f, 0.3f)) });

    std::vector<Plane> planes;
    planes.push_back( { normalize(vec3(0.f, -1.f, 1.f)), vec3(0.f), 200.f, 100.f, Material(vec3(0.1f, 0.25f, 1.f)) } );

    // ------------------- BUFFERS -------------------
    std::vector<float> depth (cam.w * cam.h);
    for(size_t i = 0; i < depth.size(); ++i)
        depth[i] = std::numeric_limits<float>::infinity();

    std::vector<unsigned char> col (cam.w * cam.h * 4);
    for(size_t i = 0; i < col.size(); i+=4){
        col[i] = 0;
        col[i+1] = 0;
        col[i+2] = 0;
        col[i+3] = 255;
    }

    // ------------------- RENDER -------------------
    for(int64_t u = 0; u < cam.w; ++u){
       for(int64_t v = 0; v < cam.h; ++v){
            vec3 ray_o, ray_d;
            
            const size_t idx_d = v * cam.w + u;
            const size_t idx_c = 4*idx_d;

            for(size_t ms = 0; ms < 3; ++ms){

                Ray ray = cameraRay(cam, u, v);

                render(col, depth, zyls, zylinder_bl, ray, ambiente, light, cam, idx_c, idx_d);
                render(col, depth, planes, plane_bl, ray, ambiente, light, cam, idx_c, idx_d);
                render(col, depth, spheres, sphere_bl, ray, ambiente, light, cam, idx_c, idx_d);

            }

        } 
    }

    std::string filename = "result_col.png";
    lodepng::encode(filename, col.data(), cam.w, cam.h);

    //std::string filename = "result_depth.png";
    //lodepng::encode(filename, col.data(), w, h);

    return EXIT_SUCCESS;
};