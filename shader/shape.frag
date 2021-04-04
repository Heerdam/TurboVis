#version 430 core

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

/*
 AND    A * B
 OR     (A + B) % 2
*/

int NOT(in int _A){
    return (_A + 1)%2;
}

int AND(in int _A, in int _B){
    return _A * _B;
}

int OR(in int _A, in int _B){
    return (_A + _B) % 2;
}

bool zylinder(in Zylinder _p, in vec3 _ray_o, in vec3 _ray_d, out vec3 _intersect, out vec3 _normal){
    const vec3 P = _p.p1;
    const vec3 Q = _p.p2;
    const float r = _p.radius;
    const vec3 n = _ray_d;
    const vec3 A = _ray_o;
    const vec3 m = A - P;
    const vec3 d = normalize(Q-P);

    const float md = dot(m, d);
    const float nd = dot(n, d);
    const float dd = dot(d, d);

/*
    if(md <= 0.f && md + nd < 0.f) return false; //e1
    if(nd > dd && dm + dn > dd) return false;  //e2
*/
    const int e1 = int(md <= 0.f && md + nd < 0.f);
    const int e2 = int(nd > dd && dm + dn > dd);

    const float nn = dot(n, n);
    const float mn = dot(m, n);

    const float a = dd*nn - nd*nd
    const float k = dot(m, m) - r*r;
    const float c = dd*k - md*md;

    float t = 0.f;

/*
    e3 = AND(e3_a, NOT(e3_b))

    if(abs(a) < EPSILON){ //e3_a
        if(c > 0.f) return false; //e3_b
        if(md < 0.f) //e3_c
            t = -mn /nn;
        else if(md > dd) //e3_d
            t = (nd - mn) /nn;
        else t = 0.f;
        return true;
    }
*/
    const int ea3_a = int(abs(a) < EPSILON);
    const int ea3_b = int(c > 0.f);
    const int ea3_c = int(md < 0.f);
    const int ea3_d = int(md > dd);

    const int e3 = AND(e3_a, NOT(e3_b));

    /*
        NOT(e3_a) -> t = 0
        AND(AND(e3_a, NOT(e3_b)), e3_c) -> t = -mn /nn
        
    */
    t = ea3_c * (-mn /nn) + (1-e3_c) * e3_d * ((nd - mn) /nn);

    const float b = dd * mn - nd * md;
    const float discr = b*b - a*c;
/*
    if(discr < 0.f) return false; //e4
*/

    const int e4 = int(discr < 0.f);
    /*
        e3 -> t

    */
    t = e3 * t + NOT(e3) * (-b - sqrt(discr)) / a;

/*
    if(t < 0.f || t > 1.f) return false; //e5
*/

    const int e5 = int(t < 0.f || t > 1.f);

/*
    //e6
    if(md + t*nd < 0.f){ //e6_a
        if(nd <= 0.f) return false; //e6_b
        t = -md/nd;
        return k + 2*t*(mm+t*nn) <= 0.f; //e6_c

    //e7
    } else if(md + t * nd > dd){ //e7_a
        if(nd >= 0.f) return false; //e7_b
        t = (dd -md)/nd;
        return k + dd - 2*md +t - (2*(mn-nd) + t*nn) <= 0.f; //e7_c
    }
*/
    const int e6_a = int(md + t*nd < 0.f);
    const int e6_b = int(nd <= 0.f);
    const int e6_c = int(k + 2*t*(mm+t*nn) <= 0.f);

    const int e7_a = int(md + t * nd > dd);
    const int e7_b = int(nd >= 0.f);
    const int e7_c =  int(k + dd - 2*md +t - (2*(mn-nd) + t*nn) <= 0.f);

    /*
        e3 -> t
        AND(AND(NOT(e3), e6_a), NOT(e6_b))  -> t = -md/nd
        AND(AND(AND(NOT(e3), NOT(e6_a)), e7_a), NOT(e7_b)) -> t = (dd -md)/nd 

    */
    t = e3 * t + (1-e3) * (e6_a * (-md/nd) + (1-e6_a) * e7_a * ((dd -md) / nd));

    const int e6 = AND(AND(e6_a, NOT(e6_b)), e6_c);
    const int e7 = AND(AND(e7_a, NOT(e7_b)), e7_c);

    //if e1 -> finito
    //if e2 -> finito
    //if e3_a -> e3
    //if e4 -> finito
    //if e5 -> finito
    //if e6_a || e7_a -> e6 || e7


    const int outv = (1-e1) * (1-e2) * ((1-e3_a) + e3) * (1-e4) * (1-e5) * ((1-e6_a) * e7_a + e6) * ((1-e7_a) * e6_a + e7);

    _intersect = _ray_o + t * _ray_d; 

    return bool(outv);

}


out vec4 fragColor;

void main() {

    const vec3 r_o;
    const vec3 r_d;

    for(int i = 0; i < sp_size; ++i){

    }

    for(int i = 0; i < zyl_size; ++i){
        
    }

    for(int i = 0; i < pln_size; ++i){
        
    }

}