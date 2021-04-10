#version 430 core

layout(local_size_x = 32​, local_size_y = 32​, local_size_z = 1​) in;

#define EPSILON 0.00001

struct Sphere {
    float radius;
    vec3 centre;
}

struct Plane {
    vec3 n;
    vec3 p0; //assumed the centre of the plane
    vec2 hs; //halfsizes
}

struct Zylinder {
    float radius;
    vec3 p1;
    vec3 p2;
}

layout(std430, binding = 1) readonly buffer sph_buf {
    int sp_size;
	Sphere spheres[];
};

layout(std430, binding = 2) readonly buffer zyl_buf {
    int zyl_size;
	Zylinder zyls[];
};

layout(std430, binding = 3) readonly buffer pln_buf {
    int pln_size;
	Plane planes[];
};

float dist2(in vec3 _v){
    return dot(_v, _v);
}

float sqr(in float _v){
    return _v*_v;
}

int NOT(in int _A){
    return (_A + 1)%2;
}

int AND(in int _A, in int _B){
    return _A * _B;
}

int OR(in int _A, in int _B){
    return (_A + _B) % 2;
}

bool sphere(in Sphere _s, in vec3 _ray_o, in vec3 _ray_d, out vec3 _intersect, out vec3 _normal) {
    const vec3 oc = _ray_o - _s.centre;
    const float d = dot(_ray_d, oc);
    const float delta = sqr(d) - (dist2(oc) - _s.radius*_s.radius);
    const float t = -dot(d) - sqrt(delta);
    _intersect = _ray_o + t*_ray_d;
    _normal = normalize(_intersect - _s.centre);
    return !(delta < 0.f);
}

bool plane(in Plane _p, in vec3 _ray_o, in vec3 _ray_d, out vec3 _intersect, out vec3 _normal){
    const float s = dot(_p.n, _p.p0 - _ray_o);
    const float t = dot(_ray_d, _p.n);
    const float d = t != 0.f ? s / t : 0.f; //todo 
    _intersect = _ray_o + d * _ray_d;
    _normal = _p.n;
    const float dist = distance(_intersect, _p.p0);
    return d != 0.f && dist <= _p.hs.x && dist <= _p.hs.y;
}

bool zylinder(in Zylinder _p, in vec3 _ray_o, in vec3 _ray_d, out vec3 _intersect, out vec3 _normal){
    const vec3 P = _p.p1;
    const vec3 Q = _p.p2;
    const float r = _p.radius;
    const vec3 n = _ray_d;
    const vec3 A = _ray_o;
    const vec3 m = A - P;
    const vec3 d = normalize(Q-P);

    const float md = Dot(m, d);
    const float nd = Dot(n, d);
    const float dd = Dot(d, d);

    // Test if segment fully outside either endcap of cylinder
    if (md < 0.f && md + nd < 0.f) return false; // Segment outside ‘p’ side of cylinder
    if (md > dd && md + nd > dd) return false;     // Segment outside ‘q’ side of cylinder

    const float nn = dot(n, n);
    const float mn = dot(m, n);
    const float a = dd * nn – nd * nd;
    const float k = dot(m, m) – r * r;
    const float c = dd * k – md * md;
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

    const float b = dd * mn – nd * md;
    const float discr = b * b – a * c;

    if (discr < 0.f) return false; // No real roots; no intersection

    t = (-b – sqrt(discr)) / a;

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
        t = (dd – md) / nd;
        // Keep intersection if Dot(S(t) - q, S(t) - q) <= r^2
        return k + dd – 2 * md + t * (2 * (mn – nd) + t * nn) <= 0.f;
    }

    // Segment intersects cylinder between the end-caps; t is correct
    return true; 
}


layout (location = 10) uniform vec3 camPosWorld;
layout (location = 11) uniform mat4 proj;
layout (location = 12) uniform vec2 bounds;
layout (location = 13) uniform float fovX;

layout (location = 20) uniform sampler2D depth;
layout (location = 21) uniform sampler2D color;

void main() {

    const float foxY = bounds.x / bounds.y * fovX;
    const float u;
    const float v;

    const float x = ((2*u - bounds.x)/bounds.x)*tan(fovX);
    const float y = ((2*v - bounds.y)/bounds.y)*tan(fovY);

    const vec3 r_o = camPosWorld;
    const vec3 r_d = proj * vec3(x, y, -1.f);

    for(int i = 0; i < zyl_size; ++i){
        vec3 intersection;
        vec3 normal;
        const bool hit = zylinder(zyls[i], r_o, r_d, intersection, normal);
        if(hit){
            const float depth = texelFetch(depth, ivec2(x, y)).z;
        }
    }


}
