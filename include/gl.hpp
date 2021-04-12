#ifndef GLL_HPP
#define GLL_HPP

#include "defines.hpp"

namespace GL {

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
	        compute shader: [PATH_TO_FILE].comp
	        vertex shader: [PATH_TO_FILE].ver
	        geometry shader: [PATH_TO_FILE].geo
	        fragment shader: [PATH_TO_FILE].frag
	    */
        bool compileFromFile(const std::string&);
        bool compile(const char*, const char*, const char*, const char*);
        GLuint getHandle() const;
        void bind() const;
        void unbind() const;
    };

    class Camera {

    public:

        const Vec3 upAxis = Vec3(0.f, 1.f, 0.f);

        Camera(int64_t _width, int64_t _height, float _fov, float _near, float _far) noexcept
        : width(_width), height(_height), fov(_fov), near(_near), far(_far) {}

        void update() noexcept{
            dir = normalize(target - position);
            right = normalize(cross(dir, upAxis));
            up = normalize(cross(right, dir));
        }

        Vec3 position, target;
        Vec3 up, right, dir;

        const float near, far, fov;
        const int64_t width, height;

        Mat4 combined;

        static Vec3 shoemake_projection(const Vec2&, float);
        static Vec3 holroyd_projection(const Vec2&, float);
        static Quat trackball_shoemake(const Vec2&, const Vec2&, float);
        static Quat trackball_holroyd(const Vec2&, const Vec2&, float);

        [[nodiscard]] Vec3& operator*(Vec3& _in) const noexcept {
            _in = normalize(_in.x*right + _in.y*up + _in.z*dir);
            return _in;
        };

    };

    class ShapeRenderer{

        struct Sphere {
            float radius;
            Vec3 centre;
            Vec3 col;
        };

        struct Plane {
            Vec3 n;
            Vec3 p0;  //assumed the centre of the plane
            float hh, hw;
            Vec3 col;
        };

        struct Zylinder {
            float radius;
            Vec3 p1;
            Vec3 p2;
            Vec3 col;
        };

        ShaderProgram shader;

        GLuint VAO; //full screen quad

        // ------------------- SSBO -------------------

        GLuint SSBO_spheres[2];
        GLuint SSBO_zyls[2];
        GLuint SSBO_planes[2];

        Sphere* SSBO_ptr_spheres[2];
        Zylinder* SSBO_ptr_zyls[2];
        Plane* SSBO_ptr_planes[2];

        size_t currIndex = 0;

        size_t spheres, zyls, planes;
        const size_t maxShapes;

    public:
        ShapeRenderer(uint32_t /*_maxShapesPerType*/);
        void render(const Camera& /*_cam*/);
        //Lines
        void drawLine(const Vec3& /*_p1*/, const Vec3& /*_p2*/, float /*_radius*/, const Vec3& /*_col*/);        
        void drawAABB(const Vec3&, const Vec3&, float, const Vec3&);
        void drawSphere(const Vec3& /*_centre*/, float /*_radius*/, const Vec3& /*_col*/);
        void drawAxisWidget();
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