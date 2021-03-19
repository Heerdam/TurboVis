#version 430 core

layout (location = 0) in vec3 pos;

layout (location = 1) uniform mat4 cam;

out vec4 col;

void main() {
    col = vec4(1., 0, 0, 0);
    gl_Position = cam * vec4(pos.xyz, 1);
}