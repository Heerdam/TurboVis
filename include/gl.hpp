#ifndef GL_HPP
#define GL_HPP

#include "defines.hpp"
#include "math.hpp"
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

    namespace util {
        struct alignas(16) array4f16a {
            float v[4];
            float operator[](size_t _index){ return v[_index]; }
        };
        struct alignas(32) array4f32a {
            float v[8];
            float operator[](size_t _index){ return v[_index]; }
        };
        struct alignas(16) array4ui16a {
            uint16_t v[16];
            uint16_t operator[](size_t _index){ return v[_index]; }
        };
        struct alignas(32) array4ui32a {
            uint32_t v[8];
            uint32_t operator[](size_t _index){ return v[_index]; }
        };
    } //namespace util

    /*

    template<class Func>
    struct MediumConfig {
        size_t subdivisions = 4;
        
        Func function; //float foo(Vec3)
        std::string function_shader;
    };

    template <class T = float, class LightDim = 26>
    class Medium {

        struct Triangle{ size_t idx[3] };

        using Vector = Eigen::Matrix<T, 3, 1>;
        using LightVec = Eigen::Matrix<T, LightDim, 1>;

        std::vector<Vector> vrt;
        std::vector<Triangle> trias;

        std::vector<std::vector<LightVec>> signal; //voxels, bases

    public:
        Medium(const MediumConfig& _config){
            icosahedron(_subdivisions);
        }

        void computeAndUpload();
        void render();

    private:
        
        void icosahedron(size_t _subdivisions) {
        
            const T X = .525731112119133606;
            const T Z = .850650808352039932;
            const T N = 0.;

            struct Triangle { IdxType vertex[3]; };

            using TriangleList = std::vector<Triangle>;
            using VertexList = std::vector<Vector>;

            const VertexList vertices =
                {
                    {-X,N,Z}, //0
                    {X,N,Z}, //1
                    {-X,N,-Z}, //2*
                    {X,N,-Z}, //3*
                    {N,Z,X}, //4
                    {N,Z,-X}, //5* -X
                    {N,-Z,X}, //6
                    {N,-Z,-X}, //7* -X
                    {Z,X,N}, //8
                    {-Z,X, N}, //9
                    {Z,-X,N}, //10
                    {-Z,-X, N} //11
                };

            const TriangleList triangles =
                {                   
                    {0,4,1},    {0,9,4},    {4,8,1},
                    {8,10,1},   {11,0,6},   {0,1,6},
                    {6,1,10},   {9,0,11},   {5,2,3},
                    {2,7,3},    {7,2,11},   {9,2,5},
                    {9,11,2},   {7,10,3},   {8,3,10},
                    {5,3,8},    {9,5,4},    {4,5,8},
                    {7,6,10},   {7,11,6},
                };

            using Lookup = std::map<std::pair<IdxType, IdxType>, IdxType>;

            const auto vertex_for_edge = [](Lookup& _lookup, VertexList& _vertices, IdxType _first, IdxType _second) {
                typename Lookup::key_type key(_first, _second);
                if (key.first > key.second)
                    std::swap(key.first, key.second);

                const auto& inserted = _lookup.insert({ key, IdxType(vertices.size()) });
                if (inserted.second)
                    _vertices.push_back(glm::normalize(vertices[_first] + vertices[_second]));

                return inserted.first->second;
            };

            const auto subdivide = [](const VertexList& _vertices, const TriangleList& _triangles) {
                Lookup lookup;
                TriangleList result;

                for (const auto& each : triangles) {
                    std::array<IdxType, 3> mid;
                    for (size_t edge = 0; edge < 3; ++edge)
                        mid[edge] = vertex_for_edge(lookup, vertices, each.vertex[edge], each.vertex[(edge + 1) % 3]);

                    result.push_back({ each.vertex[0], mid[0], mid[2] });
                    result.push_back({ each.vertex[1], mid[1], mid[0] });
                    result.push_back({ each.vertex[2], mid[2], mid[1] });
                    result.push_back({ mid[0], mid[1], mid[2] });
                }

                return result;
            };

            VertexList vert = vertices;
            TriangleList tria = triangles;

            for (size_t i = 0; i < _subdivisions; ++i) {
                tria = subdivide(vert, tria);
            }

            vrt.reserve(vert.size() * 3);
            idx.reserve(tria.size() * 3);

            for (size_t i = 0, j = 0; i < vert.size(); ++i, j += 3) {
                const auto v = glm::normalize(vert[i]);
                std::memcpy(vert.data() + j, glm::value_ptr(v), 3 * sizeof(T));
                //vrt[j] = v[0];
                //vrt[j + 1] = v[1];
                //vrt[j + 2] = v[2];
            }

            for (size_t i = 0, j = 0; i < tria.size(); ++i, j += 3) {
                const auto& t = tria[i];
                std::memcpy(idx.data() + j, glm::value_ptr(t), 3 * sizeof(IdxType));
                //idx[j] = t.vertex[0];
                //idx[j + 1] = t.vertex[1];
                //idx[j + 2] = t.vertex[2];
            }

        }

    };

    */


}  // namespace GL

#endif /* GL_HPP */
