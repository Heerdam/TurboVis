#ifndef GLL_HPP
#define GLL_HPP

#include "defines.hpp"

namespace GL {

    namespace Frustum {
        GL_Frustum extractPlanes(const GL_Matrix4&);
        bool isPointInside(const GL_Frustum&, const Eigen::Vector3f&);
    }

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

        Eigen::Vector3f position;
        Eigen::Vector3f up;
        Eigen::Vector2f size;
        float near = 0.1f, far = 1000.f;
        GL_Matrix4 combined;

        static GL_Matrix4 projection(float, float, float, float, float, float);
        static GL_Matrix4 projection(float, float, float, float);
        static GL_Matrix4 orthographic (float, float, float, float);
        static GL_Matrix4 lookAt(const Eigen::Vector3f&, const Eigen::Vector3f&, const Eigen::Vector3f&);
        static GL_Matrix4 rotate(float, const Eigen::Vector3f&);
        static GL_Matrix4 translate(const Eigen::Vector3f&);
        static Eigen::Vector3f shoemake_projection(const Eigen::Vector2f&, float);
        static Eigen::Vector3f holroyd_projection(const Eigen::Vector2f&, float);
        static Eigen::Quaternionf trackball_shoemake(const Eigen::Vector2f&, const Eigen::Vector2f&, float);
        static Eigen::Quaternionf trackball_holroyd(const Eigen::Vector2f&, const Eigen::Vector2f&, float);
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