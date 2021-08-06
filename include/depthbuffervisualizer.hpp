#ifndef DEPTHBUFFERVISUALIZER_HPP
#define DEPTHBUFFERVISUALIZER_HPP

#include <glad/glad.h>

#include<iostream>

namespace GL {

    template <class ShaderProgram, class Camera>
    class DepthBufferVisualizer {
        ShaderProgram shader;
        GLuint VAO, FB, TEX;
        const size_t w, h;
    public:
        DepthBufferVisualizer(const Camera&);
        void render();
    };

}

template <class ShaderProgram, class Camera>
GL::DepthBufferVisualizer<ShaderProgram, Camera>::DepthBufferVisualizer(const Camera& _cam) : w(_cam.width), h(_cam.height) {

    // -------------- BUFFERS --------------

    const float quad[] = {
        1.f, 1.f, 1.f, 1.f,    // top right
        1.f, -1.f, 1.f, 0.f,   // bottom right
        -1.f, -1.f, 0.f, 0.f,  // bottom left
        -1.f, 1.f, 0.f, 1.f   // top left
    };

    const int16_t quadI[] = {
        0, 1, 2,
        2, 3, 0};

    GLuint VBO, EBO;

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    //vbo
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quad), quad, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2*sizeof(float)));
    glEnableVertexAttribArray(1);

    //ebo
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(quadI), quadI, GL_STATIC_DRAW);

    glBindVertexArray(0);

    // -------------- SHADER --------------

    const char* vertex =
        "#version 430 core\n"
        "layout (location = 0) in vec2 pos;\n"
        "layout (location = 1) in vec2 uvs;\n"
        "out vec2 uv;\n"
        "void main() {\n"
            "uv = uvs;\n"
            "gl_Position = vec4(pos.xy, 0.f, 1.f);\n"
        "};\n";

    const char* frag =
        "#version 430 core\n"
        "in vec2 uv;\n"
        "out vec4 fragColor;\n"
        "layout (location = 2) uniform sampler2D depth;\n"
        "void main() {\n"
            "fragColor = vec4(texture(depth, uv).xyz, 1.f);\n"
           // "fragColor = vec4(1.f, 0.f, 0.f, 1.f);"
        "};\n";

    shader.compile("", vertex, "", frag);

    // -------------- FRAMEBUFFER --------------

    glGenTextures(1, &TEX);
    glBindTexture(GL_TEXTURE_2D, TEX);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, _cam.width, _cam.height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    glGenFramebuffers(1, &FB);
    glBindFramebuffer(GL_FRAMEBUFFER, FB);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, TEX, 0);

    GLuint status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
    if (GL_FRAMEBUFFER_COMPLETE != status)
        std::cerr << "[DepthBufferVisualizer] couldn't create depth render target" << std::endl;


    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

template <class ShaderProgram, class Camera>
GL::DepthBufferVisualizer<ShaderProgram, Camera>::render(){

    glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, FB);
    glBlitFramebuffer(0, 0, w, h, 0, 0, w, h, GL_DEPTH_BUFFER_BIT, GL_NEAREST);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    glDisable(GL_DEPTH_TEST);
    shader.bind();
    glUniform1i(2, 0);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, TEX);
    glBindVertexArray(VAO);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_SHORT, (void*)0);
    glBindVertexArray(0);
    shader.unbind();
    glEnable(GL_DEPTH_TEST);
}

#endif //DEPTHBUFFERVISUALIZER_HPP