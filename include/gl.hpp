#ifndef GL_HPP
#define GL_HPP

#include "defines.hpp"
#include "math.hpp"
#include "camera.hpp"
#include "hdf5.hpp"
#include <lodepng.h>
#include "AsyncRenderer.hpp"

namespace GL {

    class Camera;

    class ShaderProgram {
        enum class Status {
            success,
            failed,
            missing
        };

        GLint program = -1, compute = -1, vertex = -1, geom = -1, frag = -1;
        void print(std::string, Status, Status, Status, Status, Status, std::string) const;

       public:
        ShaderProgram(std::string);
        ShaderProgram();
        ~ShaderProgram();

        std::string id;
        bool printDebug = true;

        /*
	        assumes the following extensions:
	        compute shader: [FILENAME].comp
	        vertex shader: [FILENAME].ver
	        geometry shader: [FILENAME].geo
	        fragment shader: [FILENAME].frag
	    */
        bool compileFromFile(const std::string&);
        bool compile(const char*, const char*, const char*, const char*);
        GLuint getHandle() const;
        void bind() const;
        void unbind() const;
    };

    class DepthBufferVisualizer {
        ShaderProgram shader;
        GLuint VAO, FB, TEX;
        const size_t w, h;
    public:
        DepthBufferVisualizer(const Camera&);
        void render();
    };

    struct Uniforms {      
        float fovX = 0.785398f; //45°
        float tmax = 15.f;
        int32_t steps = 3000;
        Vec3 low = Vec3(0.f);
        Vec3 high = Vec3(2.f);
        float tr_fac = 1.f;
        float pd = 1.f;
        float t = 0.01f;
        bool tt = false;
        bool grayscale = true;
        bool isosurface = false;
        float isvalue = 1.f;
    };

    class RaymarchTester {
        ShaderProgram shader;
        GLuint VAO;
    public:
        RaymarchTester();
        void render(const Camera&, const Uniforms&);
    };
 
 
    class CpuRaymarcher {

        template<class T>
        [[nodiscard]] constexpr T fract(T _x) const noexcept {
            return _x - std::floor(_x);
        }

        template<class T>
        [[nodiscard]] constexpr Eigen::Matrix<T, 3, 1> toRGB(const Eigen::Matrix<T, 3, 1>& _c) const noexcept{
            using vec3 = Eigen::Matrix<T, 3, 1>;
            using vec4 = Eigen::Matrix<T, 4, 1>;
            const vec4 K = vec4(1., 2. / 3., 1. / 3., 3.);
            const vec3 p = vec3 ( 
                std::abs( fract(_c(0) + K(0)) * 6. - K(3) ),
                std::abs( fract(_c(0) + K(1)) * 6. - K(3) ),
                std::abs( fract(_c(0) + K(2)) * 6. - K(3) ) );

            return vec3(
                _c(2) * std::lerp( K(0), std::clamp( (p(0) - K(0)) , 0., 1.), _c(1) ),
                _c(2) * std::lerp( K(0), std::clamp( (p(0) - K(0)) , 0., 1.), _c(1) ),
                _c(2) * std::lerp( K(0), std::clamp( (p(0) - K(0)) , 0., 1.), _c(1) )
            );
        }

        template<class T>
        [[nodiscard]] constexpr Eigen::Matrix<T, 3, 1> toHSL(const Eigen::Matrix<T, 3, 1>& _c) const noexcept{
            T h = 0.;
            T s = 0.;
            T l = 0.;
            const T r = _c(0);
            const T g = _c(1);
            const T b = _c(2);
            const T cMin = std::min( r, std::min( g, b ) );
            const T cMax = std::max( r, std::max( g, b ) );
            l = ( cMax + cMin ) / 2.;
            if ( cMax > cMin ) {
                const T cDelta = cMax - cMin;
                s = l < .0 ? cDelta / ( cMax + cMin ) : cDelta / ( 2. - ( cMax + cMin ) );            
                if ( r == cMax )
                    h = ( g - b ) / cDelta;
                else if ( g == cMax )
                    h = 2. + ( b - r ) / cDelta;
                else 
                    h = 4. + ( r - g ) / cDelta;

                if ( h < 0.)
                    h += 6.;

                h = h / 6.;
            }
            return Eigen::Matrix<T, 3, 1>( h, s, l );
        }

        template<class T>
        [[nodiscard]] constexpr uint16_t toColorByte(T _v) const noexcept {
            return uint16_t(_v * T(255));
        }

    public:    

        template<class Func, class Cam>
        void render(const Cam& _cam, Func _func, size_t _steps, double _tmax, size_t _width, size_t _height) const noexcept {

            using Vec4 = Eigen::Matrix<double, 4, 1>;
            using Vec = Eigen::Matrix<double, 3, 1>;
            using Vec2 = Eigen::Matrix<double, 2, 1>;

            const Vec2 bounds = Vec2(_width, _height);
            std::vector<Vec4> buffer (_width * _height);

            const Vec camPos = Vec(_cam.position.x, _cam.position.y, _cam.position.z);
            const Vec dir = Vec(_cam.dir.x, _cam.dir.y, _cam.dir.z);
            const Vec right = Vec(_cam.right.x, _cam.right.y, _cam.right.z);
            const Vec up = Vec(_cam.up.x, _cam.up.y, _cam.up.z);

            spdlog::info("render start");

            std::atomic<size_t> counter;

            #pragma omp parallel 
            {
                #pragma omp for
                for(int64_t x = 0; x < _width; ++x){
                    for(int64_t y = 0; y < _height; ++y){

                        double transmission = 1.;
                        Vec col (0., 0., 0.);

                        const Vec2 uv (x, y);
                        const Vec2 pp = (2. * (uv) - bounds) / bounds(1); //eye ray
                        const Vec r_d = ( pp(0)*right + pp(1)*up + 1.5*dir ).normalized();

                        double t = 0.;
                        const double stepsize = _tmax / _steps;

                        for( size_t s = 0; s < _steps; ++s ) {
                            const Vec pos = camPos + t * r_d;
                            //const vec3 Ls_rgb = eval_wave(vec3(1.f, 1.5f, 1.f), pos, 5, 5, 5);
                            const Vec Ls_hsl = _func(pos);
                            const Vec Ls_rgb = toRGB(Ls_hsl);
                            
                            transmission *= std::exp( -Ls_hsl(2) );
                            col += transmission * Ls_rgb;
                            t += stepsize;
                        }

                        buffer[x + y* _width] = Vec4( col(0), col(1), col(2), 1. - transmission );
                        ++counter;
                        const size_t m = ((_width * _height)/100);
                        if(counter % m == 0){
                            spdlog::debug("{}%", 100./double(_width * _height) * double(counter));
                        }

                    }
                    
                }

            } 

            spdlog::info("render done");

            std::vector<unsigned char> out(buffer.size() * 4);
            for(size_t i = 0, j = 0; i < buffer.size(); ++i, j+=4){
                out[j] = toColorByte(buffer[i](0));
                out[j+1] = toColorByte(buffer[i](1));
                out[j+2] = toColorByte(buffer[i](2));
                out[j+3] = toColorByte(buffer[i](3));
            }

            std::string filename = "render.png";
            lodepng::encode(filename, out.data(), _width, _height);

            spdlog::info("Finished");

        }

    };

    template<class T>
    class HagedornRenderer {

        Eigen::Matrix<T, 3, 1> dir, up, right;
        Eigen::Matrix<T, 3, 1> camPos;

        Eigen::Matrix<T, -1, 1> lower, upper;

        ShaderProgram shader;
        GLuint VAO[2], TEX[2];
        float* tex_buffer[2];
        int64_t index = 0;
        int64_t width, height;

        AsyncRenderer renderer;

        std::vector<T> buffer;
        std::vector<float> colbuffer;

        IO::File<T> file;

        std::atomic<size_t> t;

    public:
        HagedornRenderer(const Camera& /*_cam*/) noexcept;
        void set(const IO::File<T>& /*_file*/) noexcept;
        void start() noexcept;
        void render(const Camera& /*_cam*/) noexcept;
        void stop() noexcept;
    };//HagedornRenderer

    template <class T>
    [[nodiscard]] Eigen::Matrix<T, 3, 1> c_to_HSL(const std::complex<T>& /*_c*/) noexcept;

    template <class T>
    [[nodiscard]] Eigen::Matrix<T, 3, 1> HSL_to_RGB_rad(const Eigen::Matrix<T, 3, 1>& /*_hsv*/) noexcept;

    template <class T>
    [[nodiscard]] Eigen::Matrix<T, 3, 1> HSL_to_RGB_deg(const Eigen::Matrix<T, 3, 1>& /*_hsv*/) noexcept;

}  // namespace GL


template <class T>
inline Eigen::Matrix<T, 3, 1> GL::c_to_HSL(const std::complex<T>& _c) noexcept {
    using Vector = Eigen::Matrix<T, 3, 1>;
    const T phase = std::arg(_c);
    Vector hsv(0.5 * std::fmod(phase + 2. * M_PI, 2. * M_PI) / M_PI, 1., 1.);
    const std::complex<T> modulus = std::abs(_c);

    //lightness
    hsv(2) = 2. * std::atan2(modulus.real(), 1.) / M_PI;

    //saturation
    const T l = hsv(2);
    hsv(1) = (l <= 0.5) ? 2 * l : 2. * (1. - l);
    
    return hsv;
}; //c_to_HSV

/*
    h: [0, 2pi]
    s: [0, 1]
    l: [0, 1]
    rgb: [0, 1]
*/
template <class T>
inline Eigen::Matrix<T, 3, 1> GL::HSL_to_RGB_rad(const Eigen::Matrix<T, 3, 1>& _hsv) noexcept {
    return HSL_to_RGB_deg<T>( { _hsv(0) * T( M_PI / 180.), _hsv(1), _hsv(2) } );
}; //HSV_to_RGB

/*
    h: [0, 360]
    s: [0, 1]
    l: [0, 1]
    rgb: [0, 1]
*/
template <class T>
inline Eigen::Matrix<T, 3, 1> GL::HSL_to_RGB_deg(const Eigen::Matrix<T, 3, 1>& _hsv) noexcept {

    const T H = _hsv(0);
    const T S = _hsv(1);
    const T L = _hsv(2);

    assert(0. <= H && H <= 360.);
    assert(0. <= S && S <= 1.);
    assert(0. <= L && L <= 1.);
    
    const T C = ( T(1.) - std::abs( T(2.) * L - T(1.) ) ) * S;
    const T X = C * (T(1.) - std::abs(std::fmod(H / T(60.), T(2.)) - T(1.)));
    const T m = L - C * T(0.5);

    switch(size_t(H / 60.)){
        case 0: return { C + m, X + m, m};
        case 1: return { X + m, C + m, m};
        case 2: return { m, C + m, X + m};
        case 3: return { m, X + m, C + m};
        case 4: return { X + m, m, C + m};
        case 5: return { C + m, m, X + m};
        default: return { 0., 0., 0.};
    }

}; //HSV_to_RGB

template<class T>
inline void GL::HagedornRenderer<T>::set(const IO::File<T>& _file) noexcept {
    file = _file;
};

template<class T>
inline GL::HagedornRenderer<T>::HagedornRenderer(const GL::Camera& _cam) noexcept {
    width = _cam.width;
    height = _cam.height;
    lower.resize(3);
    upper.resize(3);
    for(size_t i = 0; i < 3; ++i){
        lower(i) = -1.;
        upper(i) = 1.;
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
        "layout(std430, binding = 2) readonly buffer tex {\n"
	        "float vals[];\n"
        "};\n"
        "layout (location = 3) uniform vec2 bounds;\n"
        "void main() {\n"
            "const int idx = 4*(int(gl_FragCoord.x) + int(gl_FragCoord.y) * int(bounds.x));\n"
            "fragColor = vec4(vals[idx], vals[idx+1], vals[idx+2], vals[idx+3]);\n"
        "};\n";

    shader.compile("", vertex, "", frag);

    // -------------- BUFFER --------------
    buffer.resize(4 * _cam.width * _cam.height);
    colbuffer.resize(4 * _cam.width * _cam.height);
    //std::memset(buffer.data(), 0, buffer.size() * sizeof(T));

    // -------------- CAMERA --------------
    camPos = { _cam.position.x, _cam.position.y, _cam.position.z };
    dir = { _cam.dir.x, _cam.dir.y, _cam.dir.z };
    up = { _cam.up.x, _cam.up.y, _cam.up.z };
    right = { _cam.right.x, _cam.right.y, _cam.right.z };
};

template<class T>
inline void GL::HagedornRenderer<T>::start() noexcept{

    renderer.start(1, width, height, [&](size_t _threat_index, size_t _x, size_t _y){

        //std::cout << _index << std::endl;

        using Vector = Eigen::Matrix<T, -1, 1>;
        const size_t steps = 1000;
        Vector r_o, r_d;

        const Eigen::Matrix<T, 2, 1> bounds = Eigen::Matrix<T, 2, 1>(T(width), T(height));
        const Eigen::Matrix<T, 2, 1> uv = Eigen::Matrix<T, 2, 1>(T(_x), T(_y));
        const Eigen::Matrix<T, 2, 1> pp = (2. * uv - bounds) / bounds(1);  //eye ray

        r_o = camPos;
        r_d = (pp(0) * right + pp(1) * up + 1.5 * dir).normalized();

        T maxDist = 25.;
        const T dS = maxDist / T(steps);

        //Eigen::array<Eigen::Index, -1> _dims (file);

        Eigen::Matrix<T, 3, 1> col;
        col.setZero();

        T t = 0.;
        T transmission = 1.;
        
        //intersect bounding box for early out
        const bool hit = Math::intersect(r_o, r_d, lower, upper, maxDist, t);
        if(!hit){
            col(0) = col(1) = col(2) = 0.1;
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
                    pos(2) > upper(2))
                    return;

                //calculate basis function
                const std::vector<std::complex<T>> phis = Math::Hagedorn::compute(
                    pos,
                    0.1,
                    file.k_max,
                    file.p[0],
                    file.q[0],
                    file.Q[0],
                    file.P[0]);

                    //calculate linear combination
                    
                    std::complex<T> res (0., 0.);
                    for(size_t i = 0; i < file.Ks.size(); ++i){
                        const auto& index = file.Ks[i];

                        //flatten multi-index
                        Eigen::Index mi = index[0];
                        for(size_t d = 0; d < index.rows(); ++d){
                            mi *= file.k_max(d);
                            mi += index[d];
                        }

                        res += file.c_0[t](i) * phis[mi];

                    } 

                    //compute color
                    const auto hsv = GL::c_to_HSL(res);
                    transmission *= std::exp(hsv(2) * dS);
                    const auto rgb = GL::HSL_to_RGB_deg(hsv);
                    col += transmission * rgb;

                    if(renderer.isShutdown() || renderer.isRestart(_threat_index))
                        return;
            }
        }        

        const size_t idx = _y * width + _x;

        buffer[4 * idx] = col(0);
        buffer[4 * idx + 1] = col(1);
        buffer[4 * idx + 2] = col(2);
        buffer[4 * idx + 3] = transmission;

        colbuffer[4 * idx] = float(col(0));
        colbuffer[4 * idx + 1] = float(col(1));
        colbuffer[4 * idx + 2] = float(col(2));
        colbuffer[4 * idx + 3] = float(transmission);
    });

};

template<class T>
inline void GL::HagedornRenderer<T>::render(const GL::Camera& _cam) noexcept {
    {
        if(_cam.hasMoved){   
            camPos = { _cam.position.x, _cam.position.y, _cam.position.z };
            dir = { _cam.dir.x, _cam.dir.y, _cam.dir.z };
            up = { _cam.up.x, _cam.up.y, _cam.up.z };
            right = { _cam.right.x, _cam.right.y, _cam.right.z };
            std::lock_guard<std::mutex> lock(renderer.getMutex());
            std::memset(colbuffer.data(), 0, colbuffer.size() * sizeof(float));
            std::memset(buffer.data(), 0, buffer.size() * sizeof(T));
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

    //render tex
    glDisable(GL_DEPTH_TEST);
    shader.bind();
    glUniform2f(3, float(_cam.width), float(_cam.height));
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, TEX[index]);
    glBindVertexArray(VAO[index]);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_SHORT, (void*)0);
    glBindVertexArray(0);
    shader.unbind();
    glEnable(GL_DEPTH_TEST);

    index = (index+1)%2;
};

template<class T>
inline void GL::HagedornRenderer<T>::stop() noexcept{
    renderer.stop();
};

#endif /* GL_HPP */
