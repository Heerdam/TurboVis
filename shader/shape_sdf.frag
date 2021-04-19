#version 430 core

float pow8(in float _x) {
    float x = _x * _x; // xˆ2
    x *= x; // xˆ4
    return x * x;
}

float length8(in vec2 _v) {
    return pow(pow8(_v.x) + pow8(_v.y), 1.f / 8.f);
}

// ------------------- STRUCTS -------------------

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

struct Wheel {
    float cX, cY, cZ;
    float r, R, thickness;
    float r11, r21, r31, r12, r22, r32, r13, r23, r33;
    float colX, colY, colZ;
};

struct Hit {
    uint shape;
    uint idx;
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

layout(std430, binding = 5) readonly buffer line_buf {
	Zylinder lines[];
};

layout(std430, binding = 6) readonly buffer wheel_buf {
	Wheel wheels[];
};

// ------------------- SDF -------------------

float zylinder(in Zylinder _p, in vec3 _pos) {
    const vec3 p1 = vec3(_p.p1X, _p.p1Y, _p.p1Z);
    const vec3 p2 = vec3(_p.p2X, _p.p2Y, _p.p2Z);
    const vec3  ba = p2 - p1;
    const vec3  pa = _pos - p1;
    const float baba = dot(ba, ba);
    const float paba = dot(pa, ba);
    const float x = length(pa * baba - ba * paba) - _p.radius * baba;
    const float y = abs(paba-baba * 0.5f) - baba * 0.5f;
    const float x2 = x * x;
    const float y2 = y * y * baba;   
    const float d = (max(x, y) < 0.f) ? -min(x2, y2) : ( ((x>0.0) ? x2 : 0.0) + ((y > 0.f) ? y2 : 0.f) );   
    return sign(d) * sqrt(abs(d)) / baba;
}

float line(in Zylinder _p, in vec3 _pos) {
    const vec3 p1 = vec3(_p.p1X, _p.p1Y, _p.p1Z);
    const vec3 p2 = vec3(_p.p2X, _p.p2Y, _p.p2Z);
    const vec3  ba = p2 - p1;
    const vec3  pa = _pos - p1;
    const float h = clamp( dot(pa,ba)/dot(ba,ba), 0.f, 1.f );
    return length( pa - ba*h ) - _p.radius;
}

float sphere(in Sphere _s, in vec3 _pos) {
    const vec3 sc = vec3(_s.centreX, _s.centreY, _s.centreZ);
    return length(_pos - sc) - _s.radius;
}

float wheel(in Wheel _w, in vec3 _pos) {
    const mat3 rot = mat3(_w.r11, _w.r21, _w.r31, _w.r12, _w.r22, _w.r32, _w.r13, _w.r23, _w.r33);
    const vec3 pos = rot * _pos;
    const vec3 q = max(abs(pos) - vec3(0.f, _w.thickness, 0.f), 0.f);
    const vec4 v = vec4( max(q, 0.f), min( max(q.x, max(q.y, q.z)), 0.f ) );
    const vec3 C = vec3(_w.cX, _w.cY, _w.cZ);
    const vec2 t = vec2(length(v.xz - C.xz) - _w.R, v.y - C.y);
    return length(t) - _w.r;
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
layout (location = 20) uniform int line_size = 0;
layout (location = 21) uniform int wheel_size = 0;

//optional
layout (location = 31) uniform float fovX = 0.785398f; //45°

layout (location = 32) uniform vec3 light_dir = normalize(vec3(1.f));
layout (location = 33) uniform vec3 light_col = vec3(1.f);
layout (location = 34) uniform float light_strength = 0.75f;

layout (location = 35) uniform vec3 ambiente_col = vec3(1.f);
layout (location = 36) uniform float ambiente_strength = 0.4;

layout (location = 37) uniform uint AA = 2;
layout (location = 38) uniform float tmax = 1500.f;
layout (location = 39) uniform uint steps = 128;

// --------------------------------------

float map(in vec3 _pos, out Hit _hit) {


    float res = 1.f / 0.f;

    //spheres
    for(int i = 0; i < sphere_size; ++i){
        const float r = sphere(spheres[i], _pos);
        if(r < res){
            res = r;
            _hit.shape = 0;
            _hit.idx = i;
        }       
    }

    //lines
    for(int i = 0; i < line_size; ++i){
        const float r = line(lines[i], _pos);
        if(r < res){
            res = r;
            _hit.shape = 1;
            _hit.idx = i;
        }       
    }

    //zylinders
    for(int i = 0; i < zyl_size; ++i){
        const float r = zylinder(zyls[i], _pos);
        if(r < res){
            res = r;
            _hit.shape = 2;
            _hit.idx = i;
        }       
    }

    //wheels
    for(int i = 0; i < wheel_size; ++i){
        const float r = wheel(wheels[i], _pos);
        if(r < res){
            res = r;
            _hit.shape = 3;
            _hit.idx = i;
        }       
    }

    return res;
}

float map_n(in vec3 _pos, in Hit _hit){

    switch(_hit.shape){
        case 0: //sphere
            return sphere(spheres[_hit.idx], _pos);
        break;
        case 1: //line
            return line(lines[_hit.idx], _pos);
        break;
        case 2: //zylinder
            return zylinder(zyls[_hit.idx], _pos);
        break;
        case 3: //wheel
            return wheel(wheels[_hit.idx], _pos);
        break;
    }

    return 0.f;
}

vec3 normals(in vec3 pos, in Hit _hit) {
    // inspired by tdhooper and klems - a way to prevent the compiler from inlining map() 4 times
    vec3 n = vec3(0.f);
    #define ZERO 0
    for(int i = ZERO; i < 4; i++) {
        vec3 e = 0.5773f * (2.f * vec3((((i+3)>>1)&1),((i>>1)&1),(i&1)) - 1.f);
        n += e * map_n(pos + 0.0005f * e, _hit).x;
    }
    return normalize(n);  
}

vec3 getCol(in Hit _hit){
    switch(_hit.shape){
        case 0: //sphere
            return vec3(spheres[_hit.idx].colX, spheres[_hit.idx].colY, spheres[_hit.idx].colZ);
        break;
        case 1: //line
            return vec3(lines[_hit.idx].colX, lines[_hit.idx].colY, lines[_hit.idx].colZ);
        break;
        case 2: //zylinder
            return vec3(zyls[_hit.idx].colX, zyls[_hit.idx].colY, zyls[_hit.idx].colZ);
        break;
        case 3: //zylinder
            return vec3(wheels[_hit.idx].colX, wheels[_hit.idx].colY, wheels[_hit.idx].colZ);
        break;
    }
    return vec3(0.f);
}

// ------------------- OUT -------------------

out vec4 fragColor;

// ------------------- MAIN -------------------

void main() {

    const vec2 uv = vec2(gl_FragCoord.x, gl_FragCoord.y);

    vec3 col = vec3(0.f);

    for(uint m = 0; m < AA; m++ ){
        for(uint n = 0; n < AA; n++ ){

            const vec2 o = vec2(float(m), float(n)) / float(AA) - 0.5f; //ms
            const vec2 pp = (2.f * (uv + o) - bounds) / bounds.y; //eye ray

            //create ray
            const vec3 r_d = normalize(pp.x*right + pp.y*up + 1.5f*dir);

            //ray marching
            float t = 0.f;
            Hit hit;
            bool hasHit = false;
            for( int i = 0; i < steps; ++i ) {
                const vec3 pos = camPosWorld + t * r_d;
                const float h = map(pos, hit);
                if( h < 0.0001f || t > tmax ) {
                    hasHit = true;
                    break;
                } 
                 t += h;
            }

            //shading
            vec3 tcol = vec3(0.0);
            if( t < tmax ) {
                vec3 pos = camPosWorld + t * r_d;
                vec3 nor = normals(pos, hit);
                const float diff = max(dot(nor, light_dir), 0.f);
                const vec3 diffl = diff * light_col * light_strength;
                const vec3 Scol = getCol(hit);
                tcol = (ambiente_col * ambiente_strength + diffl) * Scol;
            }

            col += tcol;

        }
    }

    //gamma
    col /= float(AA*AA);
    const float expy = 1.f / 2.2f;
    const vec4 resg = vec4(pow(col.x, expy), pow(col.y, expy), pow(col.z, expy), 1.f);
    fragColor = resg;

}