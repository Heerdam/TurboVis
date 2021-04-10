#version 430 core

#define EPSILON 0.00001

// ------------------- HELPER FUNCTIONS -------------------

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

// ------------------- STRUCTS -------------------

struct Ray {
    vec3 orig;
    vec3 dir;
};

struct Sphere {
    float radius;
    vec3 centre;
    vec3 col;
};

struct Plane {
    vec3 n;
    vec3 p0; //assumed the centre of the plane
    float hh, hw;
    vec3 col;
};

struct Zylinder {
    float radius;
    vec3 p1;
    vec3 p2;
    vec3 col;
};

// ------------------- BUFFERS -------------------

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

// ------------------- INTERSECTORS -------------------

bool sphere(in Sphere _s, in vec3 _ray_o, in vec3 _ray_d, out vec3 _intersect, out vec3 _normal) {
   
}

bool plane(in Plane _p, in vec3 _ray_o, in vec3 _ray_d, out vec3 _intersect, out vec3 _normal){
   
}

bool zylinder(in Zylinder _p, in vec3 _ray_o, in vec3 _ray_d, out vec3 _intersect, out vec3 _normal){
    
}

// ------------------- UNIFORMS -------------------

layout (location = 10) uniform vec3 camPosWorld;
layout (location = 11) uniform mat4 proj;
layout (location = 12) uniform vec2 bounds;
layout (location = 13) uniform float fovX;

layout (location = 20) uniform sampler2D depth;
layout (location = 21) uniform sampler2D color;

// ------------------- MAIN -------------------

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
