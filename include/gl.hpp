#ifndef GLL_HPP
#define GLL_HPP

#include "defines.hpp"

namespace GL {

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

    namespace util {
        struct alignas(16) array4f16a {
            float v[4];
            float operator[](size_t _index){ return v[_index]; }
        };
        struct alignas(32) array4f32a {
            float v[8];
            float operator[](size_t _index){ return v[_index]; }
        };
    }
}

#endif /* GLL_HPP */