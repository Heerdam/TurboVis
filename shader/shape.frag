#version 430 core

float sqr(in float _v){
    return _v*_v;
}

// ------------------- STRUCTS -------------------

struct Ray {
    vec3 orig;
    vec3 dir;
};

struct Sphere {
    float radius;
    float centreX, centreY, centreZ;
    float colX, colY, colZ;
};

struct Plane {
    float nX, nY, nZ;
    float p0X, p0Y, p0Z; //assumed the centre of the plane
    float hh, hw;
    float colX, colY, colZ;
};

struct Circle {
    float nX, nY, nZ;
    float p0X, p0Y, p0Z; //assumed the centre of the plane
    float radius_outer, radius_inner;
    float colX, colY, colZ;
};

struct Zylinder {
    float radius;
    float p1X, p1Y, p1Z;
    float p2X, p2Y, p2Z;
    float colX, colY, colZ;
};

// ------------------- BUFFERS -------------------

layout(std430, binding = 1) readonly buffer sph_buf {
	Sphere spheres[];
};

layout(std430, binding = 2) readonly buffer zyl_buf {
	Zylinder zyls[];
};

layout(std430, binding = 3) readonly buffer pln_buf {
	Plane planes[];
};

layout(std430, binding = 4) readonly buffer crc_buf {
	Circle circles[];
};

layout(std430, binding = 5) readonly buffer cpzyl_buf {
	Zylinder cp_zyls[];
};

// ------------------- INTERSECTORS -------------------

int sphere(in Sphere _s, in Ray _ray, out vec3 _intersect, out vec3 _normal, out float _t) {
    const vec3 centre = vec3(_s.centreX, _s.centreY, _s.centreZ);
    const vec3 oc = _ray.orig - centre;
    const float d = dot(oc, _ray.dir);
    const float c = (dot(oc, oc) - _s.radius * _s.radius);
    const int e1 = int(c > 0.f && d > 0.f);
    const float delta = d * d - c;
    const int e2 = int(delta < 0.f);
    _t = -d - sqrt(delta);
    const int e3 = int((_t < 0.f));
    _t = (1 - e3) * _t;
    _intersect = _ray.orig + _t*_ray.dir;
    _normal = normalize(vec3(_intersect - centre));
    const int e4 = int(delta >= 0.f);
    const int res = (1 - e1)*(1 - e2)*e4;
    return res;
}

int plane(in Plane _p, in Ray _ray, out vec3 _intersect, out vec3 _normal, out float _t) {
    const vec3 n = vec3(_p.nX, _p.nY, _p.nZ);
    const vec3 p0 = vec3(_p.p0X, _p.p0Y, _p.p0Z);
    const float s = dot(n, vec3(p0 - _ray.orig));
    _t = dot(_ray.dir, n);
    const int e1 = int(_t != 0.f);
    const float d = e1 * (s / _t);
    _intersect = _ray.orig + d * _ray.dir;
    _normal = n;
    const vec3 dist = _intersect - p0;
    return int(d != 0.f && abs(dist.x) <= _p.hw && abs(dist.y) <= _p.hh);
}

int circle(in Circle _p, in Ray _ray, out vec3 _intersect, out vec3 _normal, out float _t) {
    const vec3 n = vec3(_p.nX, _p.nY, _p.nZ);
    const vec3 p0 = vec3(_p.p0X, _p.p0Y, _p.p0Z);
    const float s = dot(n, vec3(p0 - _ray.orig));
    _t = dot(_ray.dir, n);
    const int e1 = int(_t != 0.f);
    const float d = e1 * (s / _t);
    _intersect = _ray.orig + d * _ray.dir;
    _normal = n;
    const vec3 dist = _intersect - p0;
    const float dist2 = dot(dist, dist);

    return int(d != 0.f && dist2 >= sqr(_p.radius_inner) && dist2 <= sqr(_p.radius_outer));
}

int zylinder(in Zylinder _p, in Ray _ray, out vec3 _intersect, out vec3 _normal, out float _t) {
    const vec3 p1 = vec3(_p.p1X, _p.p1Y, _p.p1Z);
    const vec3 p2 = vec3(_p.p2X, _p.p2Y, _p.p2Z);
    const vec3 ba = p2 - p1;
    const vec3  oc = _ray.orig - p1;

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
    const int e2 = int(e1 == 0 && _t < 0.f);

    // body
    const float y = baoc + _t * bard;
    const int e3 = int(e1 == 0 && e2 == 0 && y > 0.f && y < baba);

    // caps !e1 && !e2 && !e3
    const int e5 = int(y < 0.f);
    const float t = ( ( (1-e5) * baba) - baoc) / bard;

    const int e4 = int(e1 == 0 && e2 == 0 && e3 == 0 && abs(k1 + k2 * t) < h);
    _t = (1-e3)*t + e3*_t;

    _intersect = _ray.orig + _t * _ray.dir;
    _normal = e3 * (normalize((oc + _t * _ray.dir - ba * y / baba) / _p.radius)) + (1-e3) * e4 * (normalize(ba * sign(y) / baba));

    const int res = int((e1 == 0 && e2 == 0) && (e3 == 1 || e4 == 1));
    return res;
}

int cp_zylinder(in Zylinder _p, in Ray _ray, out vec3 _intersect, out vec3 _normal, out float _t) {
    const vec3 p1 = vec3(_p.p1X, _p.p1Y, _p.p1Z);
    const vec3 p2 = vec3(_p.p2X, _p.p2Y, _p.p2Z);
    const vec3 ba = p2 - p1;
    const vec3  oc = _ray.orig - p1;

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
    const int e2 = int(e1 == 0 && _t < 0.f);

    // body
    const float y = baoc + _t * bard;
    const int e3 = int(e1 == 0 && e2 == 0 && y > 0.f && y < baba);

    _intersect = _ray.orig + _t * _ray.dir;
    _normal = e3 * (normalize((oc + _t * _ray.dir - ba * y / baba) / _p.radius));

    const int res = int((e1 == 0 && e2 == 0) && (e3 == 1));
    return res;
}

// ------------------- UNIFORMS -------------------

//required
layout (location = 10) uniform vec3 camPosWorld;
layout (location = 11) uniform mat4 cam;

layout (location = 12) uniform vec3 dir;
layout (location = 13) uniform vec3 right;
layout (location = 14) uniform vec3 up;

layout (location = 15) uniform vec2 bounds;

layout (location = 16) uniform int sphere_size = 0;
layout (location = 17) uniform int zyl_size = 0;
layout (location = 18) uniform int plane_size = 0;
layout (location = 19) uniform int circle_size = 0;
layout (location = 20) uniform int cp_zylinder_size = 0;

//optional
layout (location = 31) uniform float fovX = 0.785398f; //45??

layout (location = 32) uniform vec3 light_dir = normalize(vec3(1.f));
layout (location = 33) uniform vec3 light_col = vec3(1.f);
layout (location = 34) uniform float light_strength = 0.75f;

layout (location = 35) uniform vec3 ambiente_col = vec3(1.f);
layout (location = 36) uniform float ambiente_strength = 0.4;

// ------------------- OUT -------------------

out vec4 fragColor;

// ------------------- MAIN -------------------

void main() {

    // ------------------- CREATE CAMERA RAY -------------------

    const float u = gl_FragCoord.x;
    const float v = gl_FragCoord.y;

    //const float fovY = bounds.y / bounds.x * fovX;
    //const float asp = bounds.y / bounds.x;

    const float x = ((2.f*u - bounds.x)/bounds.y);// * tan(fovX);
    const float y = ((2.f*v - bounds.y)/bounds.y);// * tan(fovX);

    Ray ray;
    ray.dir = normalize(x*right + y*up + 1.5f*dir);
    ray.orig = camPosWorld;

    // ------------------- TRACING -------------------

    float t = -1.f;
    uint hit = 0;
    vec3 inter = vec3(x, y, 1.f);
    //vec3 normal;
    vec3 col = vec3(0.f);

    // ------------------- SPHERES -------------------
    for(int i = 0; i < sphere_size; ++i){
        float tt;
        vec3 tinter;
        vec3 tnor;

        const int res = sphere(spheres[i], ray, tinter, tnor, tt);
        
        const int e1 = int(res == 1 && (tt < t || t < 0.f));
        hit += e1;
        
        const uint m1 = -uint(e1 == 1);

        //t
        const uint tti = floatBitsToInt(tt);
        const uint ttir = m1 & tti;
        t = intBitsToFloat(int(ttir)) + (1 - e1) * t;

        //tinter
        const uint ttinterX = m1 & uint(floatBitsToInt(tinter.x));
        const uint ttinterY = m1 & uint(floatBitsToInt(tinter.y));
        const uint ttinterZ = m1 & uint(floatBitsToInt(tinter.z));
        inter = vec3(intBitsToFloat(int(ttinterX)), intBitsToFloat(int(ttinterY)), intBitsToFloat(int(ttinterZ))) + (1 - e1) * inter;

        const float diff = max(dot(tnor, light_dir), 0.f);
        const vec3 diffl = diff * light_col * light_strength;
        const vec3 Scol = vec3(spheres[i].colX, spheres[i].colY, spheres[i].colZ);
        const vec3 tcol = (ambiente_col * ambiente_strength + diffl) * Scol;
        
        col = e1 * tcol + (1 - e1) * col;
    }

    // ------------------- ZYLINDERS -------------------
    for(int i = 0; i < zyl_size; ++i){
        float tt;
        vec3 tinter;
        vec3 tnor;

        const int res = zylinder(zyls[i], ray, tinter, tnor, tt);

        const int e1 = int(res == 1 && (tt < t || t < 0.f));
        hit += e1;

        const uint m1 = -uint(e1 == 1);

        //t
        const uint tti = floatBitsToInt(tt);
        const uint ttir = m1 & tti;
        t = intBitsToFloat(int(ttir)) + (1 - e1) * t;

        //tinter
        const uint ttinterX = m1 & uint(floatBitsToInt(tinter.x));
        const uint ttinterY = m1 & uint(floatBitsToInt(tinter.y));
        const uint ttinterZ = m1 & uint(floatBitsToInt(tinter.z));
        inter = vec3(intBitsToFloat(int(ttinterX)), intBitsToFloat(int(ttinterY)), intBitsToFloat(int(ttinterZ))) + (1 - e1) * inter;

        const float diff = max(dot(tnor, light_dir), 0.f);
        const vec3 diffl = diff * light_col * light_strength;
        const vec3 Scol = vec3(zyls[i].colX, zyls[i].colY, zyls[i].colZ);
        const vec3 tcol = (ambiente_col * ambiente_strength + diffl) * Scol;
        
        col = e1 * tcol + (1 - e1) * col;
    }

    // ------------------- PLANES -------------------
    for(int i = 0; i < plane_size; ++i){
        float tt;
        vec3 tinter;
        vec3 tnor;

        const int res = plane(planes[i], ray, tinter, tnor, tt);

        const int e1 = int(res == 1 && (tt < t || t < 0.f));
        hit += e1;

        const uint m1 = -uint(e1 == 1);

        //t
        const uint tti = floatBitsToInt(tt);
        const uint ttir = m1 & tti;
        t = intBitsToFloat(int(ttir)) + (1 - e1) * t;

        //tinter
        const uint ttinterX = m1 & uint(floatBitsToInt(tinter.x));
        const uint ttinterY = m1 & uint(floatBitsToInt(tinter.y));
        const uint ttinterZ = m1 & uint(floatBitsToInt(tinter.z));
        inter = vec3(intBitsToFloat(int(ttinterX)), intBitsToFloat(int(ttinterY)), intBitsToFloat(int(ttinterZ))) + (1 - e1) * inter;

        const float diff = max(dot(tnor, light_dir), 0.f);
        const vec3 diffl = diff * light_col * light_strength;
        const vec3 Scol = vec3(planes[i].colX, planes[i].colY, planes[i].colZ);
        const vec3 tcol = (ambiente_col * ambiente_strength + diffl) * Scol;
        
        col = e1 * tcol + (1 - e1) * col;
    }

    // ------------------- CIRCLES -------------------
    for(int i = 0; i < circle_size; ++i){
        float tt;
        vec3 tinter;
        vec3 tnor;

        const int res = circle(circles[i], ray, tinter, tnor, tt);

        const int e1 = int(res == 1 && (tt < t || t < 0.f));
        hit += e1;

        const uint m1 = -uint(e1 == 1);

        //t
        const uint tti = floatBitsToInt(tt);
        const uint ttir = m1 & tti;
        t = intBitsToFloat(int(ttir)) + (1 - e1) * t;

        //tinter
        const uint ttinterX = m1 & uint(floatBitsToInt(tinter.x));
        const uint ttinterY = m1 & uint(floatBitsToInt(tinter.y));
        const uint ttinterZ = m1 & uint(floatBitsToInt(tinter.z));
        inter = vec3(intBitsToFloat(int(ttinterX)), intBitsToFloat(int(ttinterY)), intBitsToFloat(int(ttinterZ))) + (1 - e1) * inter;

        const float diff = max(dot(tnor, light_dir), 0.f);
        const vec3 diffl = diff * light_col * light_strength;
        const vec3 Scol = vec3(circles[i].colX, circles[i].colY, circles[i].colZ);
        const vec3 tcol = (ambiente_col * ambiente_strength + diffl) * Scol;
        
        col = e1 * tcol + (1 - e1) * col;
    }

        // ------------------- CAPELESS ZYLINDERS -------------------
    for(int i = 0; i < cp_zylinder_size; ++i){
        float tt;
        vec3 tinter;
        vec3 tnor;

        const int res = cp_zylinder(cp_zyls[i], ray, tinter, tnor, tt);

        const int e1 = int(res == 1 && (tt < t || t < 0.f));
        hit += e1;

        const uint m1 = -uint(e1 == 1);

        //t
        const uint tti = floatBitsToInt(tt);
        const uint ttir = m1 & tti;
        t = intBitsToFloat(int(ttir)) + (1 - e1) * t;

        //tinter
        const uint ttinterX = m1 & uint(floatBitsToInt(tinter.x));
        const uint ttinterY = m1 & uint(floatBitsToInt(tinter.y));
        const uint ttinterZ = m1 & uint(floatBitsToInt(tinter.z));
        inter = vec3(intBitsToFloat(int(ttinterX)), intBitsToFloat(int(ttinterY)), intBitsToFloat(int(ttinterZ))) + (1 - e1) * inter;

        const float diff = max(dot(tnor, light_dir), 0.f);
        const vec3 diffl = diff * light_col * light_strength;
        const vec3 Scol = vec3(cp_zyls[i].colX, cp_zyls[i].colY, cp_zyls[i].colZ);
        const vec3 tcol = (ambiente_col * ambiente_strength + diffl) * Scol;
        
        col = e1 * tcol + (1 - e1) * col;
    }

    const int e2 = int(hit >= 1);
    const vec4 ipos = cam * vec4(inter.xyz, 1.f);

    //calc new depth
    const float f_ndc_depth = ipos.z / ipos.w;
    const float z_ = (f_ndc_depth + 1.f) * 0.5f;
    gl_FragDepth = float(e2) * z_ + float(1-e2)*gl_FragCoord.z;

    //resulting color
    const float expy = 1.f / 2.2f;
    const vec4 resg = vec4(pow(col.x, expy), pow(col.y, expy), pow(col.z, expy), 1.f);
    fragColor = e2 * vec4(col.xyz, 1.f) + (1-e2) * vec4(0.f);
}
