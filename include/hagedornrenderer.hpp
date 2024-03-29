#ifndef HAGEDORNRENDERER_HPP
#define HAGEDORNRENDERER_HPP

#include "defines.hpp"
#include "math.hpp"
#include "camera.hpp"
#include "hdf5.hpp"
#include "shaderprogram.hpp"
#include <lodepng.h>
#include "asyncrenderer.hpp"

namespace GL {

    namespace Detail {

        

    }

    template<class T, class Camera>
    class HagedornRenderer {

        Eigen::Matrix<T, 3, 1> dir, up, right;
        Eigen::Matrix<T, 3, 1> camPos;

        Eigen::Matrix<T, -1, 1> lower, upper;

        ShaderProgram s_buffer, s_preview;
        GLuint VAO[2], TEX[2];
        GLuint pr_vao;

        float* tex_buffer[2];
        int64_t index = 0;
        int64_t width, height;

        AsyncRenderer renderer;

        std::vector<T> buffer;
        std::vector<float> colbuffer;
        //std::vector<T> depthBuffer;

        IO::File<T> file;

             

    public:
        std::atomic<size_t> steps = 500;
        std::atomic<double> k;
        std::atomic<float> scale = 5.f;
        std::atomic<float> MAX = 50.f;
        std::atomic<size_t> curT = 1;
        std::atomic<size_t> maxT = 0;
//51
    public:
        HagedornRenderer(const Camera& /*_cam*/) noexcept;
        void set(const IO::File<T>& /*_file*/) noexcept;
        void start(const Camera& /*_cam*/) noexcept;
        void render(const Camera& /*_cam*/) noexcept;
        void stop() noexcept;
        [[nodiscard]] double getProgress() noexcept;
    };//HagedornRenderer

}  // namespace GL

// --------------------- DETAIL ---------------------



// --------------------- HAGEDORNRENDERER ---------------------

template<class T, class Camera>
inline double GL::HagedornRenderer<T, Camera>::getProgress() noexcept {
    return renderer.getProgress();
}

template<class T, class Camera>
inline void GL::HagedornRenderer<T, Camera>::set(const IO::File<T>& _file) noexcept {
    file = _file;
    maxT = file.p.size();
};

template<class T, class Camera>
inline GL::HagedornRenderer<T, Camera>::HagedornRenderer(const Camera& _cam) noexcept {
    width = _cam.width;
    height = _cam.height;
    lower.resize(3);
    upper.resize(3);
    for(size_t i = 0; i < 3; ++i){
        lower(i) = -12;
        upper(i) = 12;
    }
 
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

    //vbo
    GLuint VBO;
    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quad), quad, GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    //ebo
    GLuint EBO;
    glGenBuffers(1, &EBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(quadI), quadI, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    //ssbo
    glGenBuffers(2, TEX);

    glGenVertexArrays(2, VAO);
    for(size_t i = 0; i < 2; ++i){      
        glBindVertexArray(VAO[i]);
        
        //vbo
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2*sizeof(float)));
        glEnableVertexAttribArray(1);
        
        //ebo      
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
        glBindVertexArray(0);

        //ssbo
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, TEX[i]);
        glBufferStorage(GL_SHADER_STORAGE_BUFFER, 4 * _cam.width * _cam.height * sizeof(float), nullptr, GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_CLIENT_STORAGE_BIT);
        tex_buffer[i] = reinterpret_cast<float*>(glMapBufferRange(GL_SHADER_STORAGE_BUFFER, 0, 4 * _cam.width * _cam.height * sizeof(float), GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_MAP_UNSYNCHRONIZED_BIT));
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    }

    {
        const float quad[] = {
            //x
            0.f, 0.f, 0.f, 1.f, 0.f, 0.f,   
            500.f, 0.f, 0.f, 1.f, 0.f, 0.f,  
            //y
            0.f, 0.f, 0.f, 0.f, 1.f, 0.f,   
            0.f, 500.f, 0.f, 0.f, 1.f, 0.f,
            //z
            0.f, 0.f, 0.f, 0.f, 0.f, 1.f,   
            0.f, 0.f, 500.f, 0.f, 0.f, 1.f,
 
        };

        const int16_t quadI[] = {
            0, 1, 
            2, 3,
            4, 5
        };

        GLuint pr_vbo, pr_ebo;

        glGenVertexArrays(1, &pr_vao);
        glGenBuffers(1, &pr_vbo);
        glGenBuffers(1, &pr_ebo);

        glBindVertexArray(pr_vao);

        //vbo
        glBindBuffer(GL_ARRAY_BUFFER, pr_vbo);
        glBufferData(GL_ARRAY_BUFFER, sizeof(quad), quad, GL_STATIC_DRAW);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3*sizeof(float)));
        glEnableVertexAttribArray(1);

        //ebo
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, pr_ebo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(quadI), quadI, GL_STATIC_DRAW);

        glBindVertexArray(0);

        // -------------- SHADERS --------------
        //draw buffer shader
        const char* vertex =
            "#version 430 core\n"
            "layout (location = 0) in vec2 pos;\n"
            "layout (location = 1) in vec2 uvs;\n"
            "out vec2 uv;\n"
            "void main() {\n"
            "   uv = uvs;\n"
            "   gl_Position = vec4(pos.xy, 0.f, 1.f);\n"
            "};\n";

        const char* frag =
            "#version 430 core\n"
            "in vec2 uv;\n"
            "out vec4 fragColor;\n"
            "layout(std430, binding = 2) readonly buffer tex {\n"
            "   float vals[];\n"
            "};\n"
            "layout (location = 3) uniform vec2 bounds;\n"
            "void main() {\n"
            "   const int idx = 4*( int(gl_FragCoord.x) + int(gl_FragCoord.y) * int(bounds.x) );\n"
            "   fragColor = vec4(vals[idx], vals[idx+1], vals[idx+2], vals[idx+3]);\n"
            "};\n";

        s_buffer.compile("", vertex, "", frag);
    }

    { //axis shader
        const char* vertexa = 
            "#version 430 core\n"
            "layout (location = 0) in vec3 pos;\n"
            "layout (location = 1) in vec3 col;\n"
            "layout (location = 2) uniform mat4 cam;\n"
            "out vec3 color;\n"
            "void main() {\n"
            "   color = col;\n"
            "   gl_Position =  cam * vec4(pos.xyz, 1.f);\n"
            "};\n";

        const char* fraga = 
            "#version 430 core\n"
            "in vec3 color;\n"
            "out vec4 fragColor;\n"
            "void main() {\n"
            "   fragColor = vec4(color.xyz, 1.f);\n"
             //"fragColor = vec4(1.f, 0.f, 0.f, 1.f);"
            "};\n";

        s_preview.id = "axis";
        s_preview.compile("", vertexa, "", fraga);
    }
    // -------------- BUFFER --------------
    buffer.resize(4 * _cam.width * _cam.height);
    colbuffer.resize(4 * _cam.width * _cam.height);
    //depthBuffer.resize(_cam.width * _cam.height);
    //std::memset(buffer.data(), 0, buffer.size() * sizeof(T));

    // -------------- CAMERA --------------
    camPos = { _cam.position.x, _cam.position.y, _cam.position.z };
    dir = { _cam.dir.x, _cam.dir.y, _cam.dir.z };
    up = { _cam.up.x, _cam.up.y, _cam.up.z };
    right = { _cam.right.x, _cam.right.y, _cam.right.z };
};

template<class T, class Camera>
inline void GL::HagedornRenderer<T, Camera>::start(const Camera& _cam) noexcept{

    const auto inv = Math::Hagedorn::computeInvariants(file);

    renderer.start(8, width, height, [&](size_t _threat_index, size_t _x, size_t _y){

        //std::cout << _index << std::endl;

        using Vector = Eigen::Matrix<T, -1, 1>;
        //const size_t steps = 2500;
        Vector r_o, r_d;

        const Eigen::Matrix<T, 2, 1> bounds = Eigen::Matrix<T, 2, 1>(T(width), T(height));
        const Eigen::Matrix<T, 2, 1> uv = Eigen::Matrix<T, 2, 1>(T(_x), T(_y));
        const Eigen::Matrix<T, 2, 1> pp = (2. * uv - bounds) / bounds(1);  //eye ray

        r_o = camPos + scale*(pp(0) * right + pp(1) * up );
        //r_d = (pp(0) * right + pp(1) * up + scale * dir).normalized();

        r_d = dir;

        T maxDist = 35.;
        const T dS = maxDist / T(steps);

        //Eigen::array<Eigen::Index, -1> _dims (file);

        Eigen::Matrix<T, 3, 1> col;
        col.setZero();

        T t = 0.;
        T transmission = 1.;
        
        //intersect bounding box for early out
        const bool hit = Math::intersect(r_o, r_d, lower, upper, maxDist, t); 
        t +=dS;
        if(!hit){
            
//#ifdef NDEBUG
            //col(0) = col(1) = col(2) = 0.1;
            //transmission = 0.0;
//#else
            col(0) = col(1) = col(2) = 0.0;
//#endif

        } else if(steps == 1){   

            col(0) = 1.f;
            transmission = 0.;
            
        } else {
  
            for (size_t s = 0; s < steps; ++s) {
                const Vector pos = r_o + t * r_d;
                t += dS;

                if(t > maxDist || 
                    pos(0) < lower(0) || 
                    pos(1) < lower(1) ||
                    pos(2) < lower(2) ||
                    pos(0) > upper(0) ||
                    pos(1) > upper(1) ||
                    pos(2) > upper(2)
                )
                    break;

                if constexpr (true){
                    //calculate basis function
                    
                    const std::unordered_map<Eigen::Index, std::complex<double>> phis = Math::Hagedorn::compute(
                        curT,
                        pos,
                        inv
                    );
                        

                    //calculate linear combination 
                    /*    
                    std::complex<T> res = Math::Hagedorn::Detail::phi_0(
                        pos,
                        1.,
                        file.p[curT],
                        file.q[curT],
                        file.Q[curT],
                        file.P[curT]
                    );
                    */

                    std::complex<T> res (0., 0.);
                    for(Eigen::Index k = 0; k < file.Ks.size(); ++k){ 
                        const Eigen::Index idx = Math::Hagedorn::Detail::index(file.Ks[k], file.k_max);
                        assert(phis.contains(idx));
                        const std::complex<double>& p = (*phis.find(idx)).second;
                        res += file.c_0[curT](k) * p;
                    } 
                    res *= std::exp(std::complex<T>(0., 1.) * file.S(Eigen::Index(curT)));

                    //if(res.real() > std::numeric_limits<T>::epsilon() || res.imag() > std::numeric_limits<T>::epsilon())
                       //std::cout << res << std::endl;

                    //compute color
                    const auto hsl = Math::c_to_HSL(T(MAX), res);
                    transmission *= std::exp(-hsl(2) * dS);
                    const auto rgb = Math::HSL_to_RGB_rad(hsl);
                    col += transmission * rgb * dS;

                } else {

                    Eigen::Matrix<T, 3, 1> rgb;
                    rgb.setZero();
                    rgb(0) = rgb(1) = rgb(2) = std::clamp(std::abs(std::cos(pos(0) * pos(1) * pos(2)*250. )), 0., 1.);
                    transmission *= std::exp(-rgb(2) * dS);
                    col += transmission * rgb * dS;

                }

                if(renderer.isShutdown() || renderer.isRestart(_threat_index))
                    return;
            }
        }        

        const size_t idx = _y * width + _x;

        //if(hit)
            //std::cout << col(0) << ", " << col(1) << ", " << col(2) <<  ", " << transmission << std::endl;

        buffer[4 * idx] = col(0);
        buffer[4 * idx + 1] = col(1);
        buffer[4 * idx + 2] = col(2);
        buffer[4 * idx + 3] = 1. - transmission;

        colbuffer[4 * idx] = float(col(0));
        colbuffer[4 * idx + 1] = float(col(1));
        colbuffer[4 * idx + 2] = float(col(2));
        colbuffer[4 * idx + 3] = 1.f - float(transmission);
    });

};

template<class T, class Camera>
inline void GL::HagedornRenderer<T, Camera>::render(const Camera& _cam) noexcept {
    if(true){
        if(_cam.hasMoved){   
            camPos = { _cam.position.x, _cam.position.y, _cam.position.z };
            dir = { _cam.dir.x, _cam.dir.y, _cam.dir.z };
            up = { _cam.up.x, _cam.up.y, _cam.up.z };
            right = { _cam.right.x, _cam.right.y, _cam.right.z };
            //std::lock_guard<std::mutex> lock(renderer.getMutex());
            std::memset(colbuffer.data(), 0, colbuffer.size() * sizeof(float));
            std::memset(buffer.data(), 0, buffer.size() * sizeof(T));
            //for(size_t i = 0; i < depthBuffer; ++i)
                //depthBuffer[i] = std::numeric_limits<T>::infinity();
            renderer.restart();
        }

        //upload to buffer
        std::memcpy(tex_buffer[index], colbuffer.data(), colbuffer.size() * sizeof(float));
        /*
        for(size_t i = 0; i < width * height; ++i){
            tex_buffer[index][4*i] = 1.f;
            tex_buffer[index][4*i+1] = 1.f;
            tex_buffer[index][4*i+2] = 0.f;
            tex_buffer[index][4*i+3] = 1.f;
        }
        */

    }

    //render axes
    //glDisable(GL_DEPTH_TEST);
    glLineWidth(5.f);
    s_preview.bind();
    glBindVertexArray(pr_vao);
    glUniformMatrix4fv(2, 1, false, glm::value_ptr(_cam.comb));
    glDrawElements(GL_LINES, 6, GL_UNSIGNED_SHORT, (void*)0);
    glBindVertexArray(0);
    s_preview.unbind();
    

    //render tex
    
    glDisable(GL_DEPTH_TEST);
    s_buffer.bind();
    glUniform2f(3, float(_cam.width), float(_cam.height));
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, TEX[index]);
    glBindVertexArray(VAO[index]);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_SHORT, (void*)0);
    glBindVertexArray(0);
    s_buffer.unbind();
    glEnable(GL_DEPTH_TEST);

    index = (index+1)%2;
};

template<class T, class Camera>
inline void GL::HagedornRenderer<T, Camera>::stop() noexcept{
    renderer.stop();
};

#endif /* HAGEDORNRENDERER_HPP */
