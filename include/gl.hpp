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

    struct Camera {

        Camera();

        Vec3 position;
        Vec3 up;
        Vec2 size;
        float near = 0.1f, far = 1000.f;
        Mat4 combined;

        static Vec3 shoemake_projection(const Vec2&, float);
        static Vec3 holroyd_projection(const Vec2&, float);
        static Quat trackball_shoemake(const Vec2&, const Vec2&, float);
        static Quat trackball_holroyd(const Vec2&, const Vec2&, float);
    };

    class ShapeRenderer{

        ShaderProgram shader;

        GLuint VAO[2], EBO[2];
        GLuint VBO_POS[2], VBO_COL[2], VBO_NRM[2];

        float* VBO_ptr_pos[2];
        float* VBO_ptr_col[2];
        float* VBO_ptr_nrm[2];
        uint16_t* EBO_ptr[2];

        uint32_t currIndex = 0;

        uint32_t currVert = 0;
        uint32_t currIdx = 0;

        std::vector<std::optional<std::pair<Vector_af32, Vector_aui16>>> spheres;

    public:
        ShapeRenderer(uint32_t /*_maxVertices*/);
        void render(const float* /*_camera*/);
        //Lines
        void drawLine(const Vec3& /*_p1*/, const Vec3& /*_p2*/, uint32_t /*_segments*/, float /*_thickness*/, const Vec4& /*_col*/);        
        void drawAABB(const Vec3&, const Vec3&, float, const Vec4&);
        void drawSphere(const Vec3& /*_centre*/, float /*_radius*/, uint32_t /*_subdivisions*/, const Vec4& /*_col*/);
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
        struct Geometry {
	        static std::pair<Vector_af32, Vector_aui16>Icosahedron(uint16_t /*_subdivisions*/);
            static std::pair<Vector_af32, Vector_aui16>Zylinder(uint16_t /*_subdivisions*/);
        };
    }
}

#endif /* GLL_HPP */