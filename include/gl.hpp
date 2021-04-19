#ifndef GLL_HPP
#define GLL_HPP

#include "defines.hpp"
#include "math.hpp"

namespace GL {

    class Camera;

    using Vector_af32 = std::vector<__m256, aligned_allocator<__m256, sizeof(__m256)>>;
    using Vector_aui16 = std::vector<__m256i, aligned_allocator<__m256i, sizeof(__m256i)>>;

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
        void drawAxisWidget() noexcept;
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



}  // namespace GL

#endif /* GLL_HPP */
