#version 430 core

#define M_PI 3.141592
#define EPSILON 0.000001

vec3 toRGB(in vec3 c){
    vec4 K = vec4(1., 2. / 3., 1. / 3., 3.);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6. - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0., 1.), c.y);
};

vec3 toHSL(in vec3 c){
    float h = 0.;
    float s = 0.;
    float l = 0.;
    float r = c.r;
    float g = c.g;
    float b = c.b;
    float cMin = min( r, min( g, b ) );
    float cMax = max( r, max( g, b ) );

    l = ( cMax + cMin ) / 2.0;
    if ( cMax > cMin ) {
        float cDelta = cMax - cMin;
        
        //s = l < .05 ? cDelta / ( cMax + cMin ) : cDelta / ( 2.0 - ( cMax + cMin ) ); Original
        s = l < .0 ? cDelta / ( cMax + cMin ) : cDelta / ( 2.0 - ( cMax + cMin ) );
        
        if ( r == cMax ) {
            h = ( g - b ) / cDelta;
        } else if ( g == cMax ) {
            h = 2.0 + ( b - r ) / cDelta;
        } else {
            h = 4.0 + ( r - g ) / cDelta;
        }

        if ( h < 0.0) {
            h += 6.0;
        }
        h = h / 6.0;
    }
    return vec3( h, s, l );
};

float henyey_greenstein(in float _theta, in float _g){
    const float a = 1.f / (4.f * M_PI);
    const float g2 = _g * _g;
    const float b = 1.f - g2;
    const float c = 1.f + g2 - 2.f * cos(_theta);
    return a * b / pow(c, 1.5f);
};

bool IntersectRayAABB(in vec3 _o, in vec3 _dir, in vec3 _low, in vec3 _high, in float _tmax, out float _t) {
    _t = -1.f / 0.f;

    for (int i = 0; i < 3; ++i) {
        if (abs(_dir[i]) < EPSILON) {
            if (_o[i] < _low[i] || _o[i] > _high[i]) return false;
        } else {

            const float ood = 1.f / _dir[i];
            float t1 = (_low[i] - _o[i]) * ood;
            float t2 = (_high[i] - _o[i]) * ood;

            if (t1 > t2){
                const float t = t2;
                t2 = t1;
                t1 = t;
            }

            _t = max(_t, t1);
            float tmax = min(_tmax, t2);

            if (_t > _tmax) return false;
        }
    }

    return true;
};

// ------------------- UNIFORMS -------------------
//required
layout (location = 10) uniform vec3 camPosWorld;
layout (location = 11) uniform mat4 cam;

layout (location = 12) uniform vec3 dir;
layout (location = 13) uniform vec3 right;
layout (location = 14) uniform vec3 up;

layout (location = 15) uniform vec2 bounds;

//optional
layout(location = 16) uniform vec3 low = vec3(0.f);
layout(location = 17) uniform vec3 high = vec3(2.f);

layout (location = 31) uniform float fovX = 0.785398f; //45°

layout (location = 32) uniform vec3 light_dir = normalize(vec3(1.f));
layout (location = 33) uniform vec3 light_col = vec3(1.f);
layout (location = 34) uniform float light_strength = 0.75f;

layout (location = 35) uniform vec3 ambiente_col = vec3(1.f);
layout (location = 36) uniform float ambiente_strength = 0.4;

layout (location = 37) uniform uint AA = 1;
layout (location = 38) uniform float tmax = 15.f;
layout (location = 39) uniform uint steps = 2000;

layout (location = 50) uniform float t = 0.f;
// ------------------- OUT -------------------
out vec4 fragColor;

// ------------------- FUNCTION -------------------
vec3 eval_exp(in vec3 _pos){
    if(abs(_pos.x) > 1.f || abs(_pos.y) > 1.f || abs(_pos.z) > 1.f)
        return vec3(0.f);
    const float res = clamp((exp((_pos.x) * (_pos.y) * (_pos.z))-1.f) / 3.f, -1.f, 1.f);
    if(res >= 0)
        return vec3(0, res, 0);
    else
        return vec3(abs(res), 0, 0);
};

vec3 eval_sin(in vec3 _pos){
    if(abs(_pos.x) > 1.f || abs(_pos.y) > 1.f || abs(_pos.z) > 1.f)
        return vec3(0.f);
    const float res = sin(_pos.x *t) * sin(_pos.y * t) * sin( _pos.z * t);
    if(res >= 0)
        return vec3(0, res, 0);
    else
        return vec3(abs(res), 0, 0);  
};

vec3 eval_wave(in vec3 _bounds, in vec3 _pos, in int _n1, int _n2, int _n3){
    if(abs(_pos.x) > abs(_bounds.x) || abs(_pos.y) > abs(_bounds.y) || abs(_pos.z) > abs(_bounds.z)) 
        return vec3(0.f);

    const float x = sin(_n1 * M_PI / _bounds.x * _pos.x);
    const float y = sin(_n2 * M_PI / _bounds.y * _pos.y);
    const float z = sin(_n3 * M_PI / _bounds.z * _pos.z);

    const float res = sqrt(8.f / (_bounds.x * _bounds.y * _bounds.z)) * x * y * z;
    if(res >= 0)
        return vec3(0, res, 0);
    else
        return vec3(abs(res), 0, 0);
};

vec3 eval_c(in vec3 _pos){
    const float x = _pos.x * 2.f;
    const float y = _pos.y * 2.f;
    const float z = _pos.z * 2.f;   

    const vec3 p = vec3(1.f, 1.f, 1.f);
    const vec3 k = p;
    const float kx = dot(k, _pos);
    const float k2 = dot(k, k);

    const float r = 1.;
    const float phi = mod((kx - k2 / 2.f * t)/2*M_PI, 2*M_PI);
    const float R = 0.;
    const float theta = M_PI - 2. * atan(r / R);

    const float H = phi;
    const float L = 1. - theta/M_PI;
    const float S = L;

    return vec3(H, S, L); 
};

void main() {

    const vec2 uv = vec2(gl_FragCoord.x, gl_FragCoord.y);

    vec3 col = vec3(0.f); //Hue saturation lightness
    float transmission = 1.f;

    const vec2 pp = (2.f * (uv) - bounds) / bounds.y; //eye ray

    //create ray
    const vec3 r_d = normalize(pp.x*right + pp.y*up + 1.5f*dir);

    //intersect aabb
    float t = 0.f;
    const bool hit = IntersectRayAABB(camPosWorld, r_d, low, high, tmax, t);

    if(!hit)
        return;

    t -= EPSILON;

    //ray marching
    const float stepsize = tmax / steps;

    for( int i = 0; i < steps; ++i ) {
        const vec3 pos = camPosWorld + t * r_d;
        const vec4 cpos = cam * vec4(pos, 1.f);
        const float zd = cpos.z / cpos.w;
        const float z = (zd + 1.f) * 0.5f;      

        if(/*z > gl_FragCoord.z && */pos.x <= high.x && pos.y <= high.y && pos.z <= high.z && pos.x >= low.x && pos.y >= low.y && pos.z >= low.z){
            const vec3 Ls_hsl = eval_c(pos);
            const vec3 Ls_rgb = toRGB(Ls_hsl);
            transmission *= exp(-Ls_hsl.z);
            col += transmission * Ls_rgb;
        }

        t += stepsize;
    }

    //gamma
    col /= float(AA*AA);
    const float expy = 1.f / 2.2f;
    const vec4 resg = vec4(pow(col.x, expy), pow(col.y, expy), pow(col.z, expy), 1.f - transmission);
    fragColor = resg;

}