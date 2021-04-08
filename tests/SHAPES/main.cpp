
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

struct Sphere {
    float radius;
    vec3 centre;
};

struct Plane {
    vec3 n;
    vec3 p0; //assumed the centre of the plane
    vec2 hs; //halfsizes
};

struct Zylinder {
    float radius;
    vec3 p1;
    vec3 p2;
};

bool sphere(const Sphere& _s, const vec3& _ray_o, const vec3& _ray_d, vec3& _intersect, vec3& _normal) {
    const vec3 oc = _ray_o - _s.centre;
    const float d = -dot(_ray_d, oc);
    const float delta = sqr(d) - (dist2(oc) - sqr(_s.radius));
    const float t = d - sqrt(delta);
    _intersect = _ray_o + t*_ray_d;
    _normal = normalize(vec3(_intersect - _s.centre));
    return delta >= 0.f;
};

bool plane(const Plane& _p, const vec3& _ray_o, const vec3& _ray_d, vec3& _intersect, vec3& _normal){
    const float s = dot(_p.n, vec3(_p.p0 - _ray_o));
    const float t = dot(_ray_d, _p.n);
    const float d = t != 0.f ? s / t : 0.f; //todo 
    _intersect = _ray_o + d * _ray_d;
    _normal = _p.n;
    const float dist = distance(_intersect, _p.p0);
    return d != 0.f && dist <= _p.hs[0] && dist <= _p.hs[1];
};

bool zylinder(const Zylinder& _p, const vec3& _ray_o, const vec3& _ray_d, vec3& _intersect, vec3& _normal){
    const vec3 P = _p.p1;
    const vec3 Q = _p.p2;
    const float r = _p.radius;
    const vec3 n = _ray_d;
    const vec3 A = _ray_o;
    const vec3 m = A - P;
    const vec3 d = normalize(vec3(Q-P));

    const float md = dot(m, d);
    const float nd = dot(n, d);
    const float dd = dot(d, d);

    // Test if segment fully outside either endcap of cylinder
    if (md < 0.f && md + nd < 0.f) return false; // Segment outside ‘p’ side of cylinder
    if (md > dd && md + nd > dd) return false;     // Segment outside ‘q’ side of cylinder

    const float nn = dot(n, n);
    const float mn = dot(m, n);
    const float a = dd * nn - nd * nd;
    const float k = dot(m, m) - r * r;
    const float c = dd * k - md * md;
    float t = 0.f;

    if (abs(a) < EPSILON) {
        // Segment runs parallel to cylinder axis
        if (c > 0.f) return false; // ‘a’ and thus the segment lie outside cylinder
        // Now known that segment intersects cylinder; figure out how it intersects
        if (md < 0.f) t = -mn / nn; // Intersect segment against ‘p’ endcap
        else if (md > dd) t = (nd - mn) / nn; // Intersect segment against ‘q’ endcap
        else t = 0.f; // ‘a’ lies inside cylinder
        return true;
    }

    const float b = dd * mn - nd * md;
    const float discr = b * b - a * c;

    if (discr < 0.f) return false; // No real roots; no intersection

    t = (-b - sqrt(discr)) / a;

    if (t < 0.f || t > 1.f) return false; // Intersection lies outside segment

    if (md + t * nd < 0.f) {
        // Intersection outside cylinder on ‘p’ side
        if (nd <= 0.f) return false; // Segment pointing away from endcap
        t = -md / nd;
        // Keep intersection if Dot(S(t) - p, S(t) - p) <= r^2
        return k + 2 * t * (mn + t * nn) <= 0.f;
    } else if (md + t * nd > dd) {
        // Intersection outside cylinder on ‘q’ side
        if (nd >= 0.f) return false; // Segment pointing away from endcap
        t = (dd - md) / nd;
        // Keep intersection if Dot(S(t) - q, S(t) - q) <= r^2
        return k + dd - 2 * md + t * (2 * (mn - nd) + t * nn) <= 0.f;
    }

    // Segment intersects cylinder between the end-caps; t is correct
    return true; 
};

void cameraRay(vec3& _ray_o, vec3& _ray_d, int64_t _w, int64_t _h, int64_t _u, int64_t _v, float _fovX, float _fovY, const mat4& _cam, const vec3& _camPos) {
    const float x = (static_cast<float>(2*_u - _w)/static_cast<float>(_w))*std::tan(_fovX);
    const float y = (static_cast<float>(2*_v - _h)/static_cast<float>(_h))*std::tan(_fovY);
    _ray_o = _camPos;
    const vec4 t = normalize(_cam * vec4(x, y, -1.f, -1.f));
    _ray_d = vec3(t.x, t.y, t.z);
};

vec3 mul(const mat4& _cam, const vec3& _vec){
    const vec4 t = _cam * vec4(_vec.x, _vec.y, _vec.z, 1.f);
    return vec3(t.x, t.y, t.z);
}

int main() {

    const int64_t w = 1000;
    const int64_t h = 1000;

    const float fovX = M_PI / 4;
    const float fovY = w / h * fovX;

    const vec3 camPos = vec3(0.f, 0.f, -500.f);
    const mat4 cam = glm::lookAt(camPos, vec3(0.f, 0.f, 0.f), vec3(0.f, 0.1f, 0.f));
    const mat4 cami = inverse(cam);
    
    std::vector<Zylinder> zyls;
    zyls.push_back({500.f, vec3(0.f, -250.f, 0.f), vec3(0.f, 250.f, 0.f)});

    std::vector<Sphere> spheres;
    spheres.push_back({1.3f, vec3(0.f)});

    std::vector<Plane> planes;

    std::vector<float> depth (w * h);
    std::vector<unsigned char> col (w * h * 4);
    std::memset(col.data(), 200, col.size());

    const vec3 lightDir = normalize(vec3(-1.f, -1.f, -1.f));
    float ambientStrength = 0.4;
    const vec3 ambient = ambientStrength * vec3(0.f, 0.5f, 1.f);
    
    for(int64_t u = 0; u < w; ++u){
       for(int64_t v = 0; v < h; ++v){
            vec3 ray_o, ray_d;
            cameraRay(ray_o, ray_d, w, h, u, v, fovX, fovY, cami, camPos);

            const size_t idx_d = v * w + u;
            const size_t idx_c = 4*idx_d;

            for(size_t k = 0; k < spheres.size(); ++k){
                vec3 inter;
                vec3 normal;
                if(sphere(spheres[k], ray_o, ray_d, inter, normal)){
                    const float diff = std::max(dot(normal, lightDir), 0.f);
                    const vec3 diffl = diff * vec3(1.f);
                    const vec3 res = (ambient + diffl) * vec3(65.f / 255.f, 178.f / 255.f , 210.f / 255.f);
                    const vec3 resg = glm::pow(res, vec3(1.f/2.2f));
    
                    //std::cout << idx_c << std::endl;
                    col[idx_c] = (unsigned char)(std::clamp(resg.x * 255.f, 0.f, 255.f));
                    col[idx_c + 1] = (unsigned char)(std::clamp(resg.y * 255.f, 0.f, 255.f));
                    col[idx_c + 2] = (unsigned char)(std::clamp(resg.z * 255.f, 0.f, 255.f));

                    col[idx_c + 3] = 255;
                } else {
                    col[idx_c] = 255;
                    col[idx_c + 1] = 255;
                    col[idx_c + 2] = 255;
                    col[idx_c + 3] = 255;
                }
            }

            

        } 
    }
    
    std::string filename = "test.png";
    lodepng::encode(filename, col.data(), w, h);

    return EXIT_SUCCESS;
};