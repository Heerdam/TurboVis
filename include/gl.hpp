#ifndef GL_HPP
#define GL_HPP

#include "defines.hpp"
#include "math.hpp"
#include "hdf5.hpp"
#include <lodepng.h>

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

    class ShapeRenderer {

        ShaderProgram shader;

        GLuint VAO; //full screen quad

        // ------------------- SSBO -------------------

        GLuint SSBO_spheres[2];
        GLuint SSBO_zyls[2];
        GLuint SSBO_planes[2];
        GLuint SSBO_circles[2];
        GLuint SSBO_lines[2];
        GLuint SSBO_wheels[2];

        Math::Intersector::Sphere<float>* SSBO_ptr_spheres[2];
        Math::Intersector::Zylinder<float>* SSBO_ptr_zyls[2];
        Math::Intersector::Plane<float>* SSBO_ptr_planes[2];
        Math::Intersector::Circle<float>* SSBO_ptr_circles[2];
        Math::Intersector::Zylinder<float>* SSBO_ptr_lines[2];
        Math::Intersector::Wheel<float>* SSBO_ptr_wheels[2];

        size_t currIndex = 0;

        size_t spheres, zyls, planes, circles, lines, wheels;
        const size_t maxShapes;

    public:
        ShapeRenderer(uint32_t /*_maxShapesPerType*/) noexcept;
        void render(const Camera& /*_cam*/) noexcept;
        //draw
        void drawLine(const Vec3& /*_p1*/, const Vec3& /*_p2*/, float /*_radius*/, const Vec3& /*_col*/) noexcept;  
        void drawZylinder(const Vec3& /*_p1*/, const Vec3& /*_p2*/, float /*_radius*/, const Vec3& /*_col*/) noexcept; 
        void drawSphere(const Vec3& /*_centre*/, float /*_radius*/, const Vec3& /*_col*/) noexcept;
        void drawWheel(const Vec3& /*_centre*/, float /*_r*/, float /*_R*/, float /*_thickness*/, const Vec3& /*_col*/, const Mat3& /*_rot*/) noexcept;
        
        void drawAxisWidget(const Mat3&) noexcept;

        //
        [[nodiscard]] bool cameraRayIntersect(const Math::Intersector::Ray<float>&) const noexcept;
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

        ShaderProgram shader;
        GLuint VAO[2], TEX[2];
        uint16_t* tex_buffer[2];
        size_t index = 0;
        size_t width, height;

        std::future<void> work;
        std::mutex mutex;
        bool isReset = false;
        bool isTerminate = false;

        std::vector<Eigen::Matrix<T, 4, 1>> buffer;
        std::vector<unsigned char> colbuffer;

        IO::File<T> file;

        std::atomic<size_t> t;

    public:
        HagedornRenderer(const Camera& /*_cam*/) noexcept;
        void set(const IO::File<T>& /*_file*/) noexcept;
        void start() noexcept;
        void render(const Camera& /*_cam*/) noexcept;
        void stop() noexcept;
    };//HagedornRenderer

}  // namespace GL

template<class T>
inline void GL::HagedornRenderer<T>::set(const IO::File<T>& _file) noexcept{
    file = _file;
};

template<class T>
inline GL::HagedornRenderer<T>::HagedornRenderer(const GL::Camera& _cam) noexcept {
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
        glBufferStorage(GL_SHADER_STORAGE_BUFFER, 4 * _cam.width * _cam.height, nullptr, GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_CLIENT_STORAGE_BIT);
        tex_buffer[i] = reinterpret_cast<uint16_t*>(glMapBufferRange(GL_SHADER_STORAGE_BUFFER, 0, 4 * _cam.width * _cam.height, GL_MAP_WRITE_BIT | GL_MAP_PERSISTENT_BIT | GL_MAP_UNSYNCHRONIZED_BIT));
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
        "layout (location = 2) uniform sampler2D tex;\n"
        "void main() {\n"
            "fragColor = vec4(texture(tex, uv).xyz, 1.f);\n"
        "};\n";

    shader.compile("", vertex, "", frag);

    // -------------- BUFFER --------------
    buffer.resize(4 * _cam.width * _cam.height);
    std::memset(buffer.data(), 0, buffer.size());
};

template<class T>
void inline GL::HagedornRenderer<T>::start() noexcept{
    work = std::async(std::launch::async, [&](){
        #pragma omp parallel 
        {
            #pragma omp for collapse(2)
            for(size_t x = 0; x < width; ++x){
                for(size_t y = 0; y < height; ++y){

                    using Vector = Eigen::Matrix<T, -1, 1>;

                    const size_t steps = 2000;
                    Vector r_o, r_d;
                    T maxDist = 10.;
                    const T dS = maxDist / T(steps);

                    Eigen::array<Eigen::Index, -1> _dims (file)

                    Eigen::Matrix<T, 3, 1> col;
                    col.setZero();

                    T transmission = 0.;
                    for(size_t s = 0; s < steps; ++s){

                        Vector pos = r_o + steps * dS * r_d;

                        //calculate basis function
                        const auto phis = Math::Hagedorn::compute(
                            pos,
                            0.1,
                            file.k_max,
                            file.p,
                            file.q,
                            file.Q,
                            file.P
                        );

                        //calculate linear combination
                        std::complex<T> res (0., 0.);
                        for(size_t i = 0; i < file.Ks.size(); ++i){
                            const auto& index = file.Ks(i);

                            //flatten multi-index
                            Eigen::Index mi = index[0];
                            for(size_t d = 0; d < index.rows(); ++d){
                                mi *= file.k_max(d);
                                mi += index(d);
                            }

                            res += file.c_0[t](i) * phis(mi);

                        }

                        //compute color
                        const auto res = Math::c_to_HSV(res);
                        transmission *= std::exp(res(2) * dS);
                        const auto rgb = Math::HSV_to_RGB(res);
                        col += transmission * rgb;

                    }
               
                    std::lock_guard<std::mutex> lock(mutex);
                    buffer[x*y](0) *= col[0];
                    buffer[x*y](1) *= col[1];
                    buffer[x*y](2) *= col[2];
                    buffer[x*y](3) = transmission;
                }
            }
        }
    });
};

template<class T>
void inline GL::HagedornRenderer<T>::render(const GL::Camera& _cam) noexcept{
    {
        std::lock_guard<std::mutex> lock(mutex);
        isReset = _cam.hasMoved;

        if(isReset)
            std::memset(buffer.data(), 0, buffer.size());

        //upload to buffer
        std::memcpy(tex_buffer[index], buffer.data(), buffer.size());
    }

    //render tex
    glDisable(GL_DEPTH_TEST);
    shader.bind();
    glUniform1i(2, 0);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, TEX[index]);
    glBindVertexArray(VAO[index]);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_SHORT, (void*)0);
    glBindVertexArray(0);
    shader.unbind();
    glEnable(GL_DEPTH_TEST);

    index = (index+1)%2;
};

template<class T>
void inline GL::HagedornRenderer<T>::stop() noexcept{

};

#endif /* GL_HPP */
