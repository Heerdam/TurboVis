#version 430 core

#define M_PI 3.141592
#define EPSILON 0.000001

// ----------------- STRUCTS -----------------

//a complex number
struct C {
    float x, y;
};

// ----------------- MATRIX -----------------

//complex conjugated of a matrix
C[9] M_conj(in C[9] _matrix){

};

//transpose with complexe conjugation
C[9] M_t_C(in C[9] _matrix){

};

//transpose without complexe conjugation
C[9] M_t(in C[9] _matrix){

};

//multiplication of 2 complex matrices
C[9] M_mul(in C[9] _m1, in C[9] _m2){

};

//matrix exponential of an complex matrix
C[9] M_exp(in C[9] _m, in float _e){

};

//determinant of a complex matrix
C M_det(in C[9] _m){

};

//the inverse of a complex matrix
C M_inv(in C[9] _m){

};

// ----------------- VECTOR -----------------

//dot product of 2 complex vectors
C V_dot(in C[3] _v1, in C[3] _v2){

};

// ----------------- VECTOR -----------------

C wave(in vec3 _pos, in vec3 _q, in vec3 _p, in C[9] _Q, in C[9] _P){

};

// ----------------- UNIFORMS -----------------

void main() {


}