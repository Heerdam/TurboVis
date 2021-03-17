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
        void print(std::string, Status, Status, Status, Status, Status, std::string);

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
        GLuint getHandle();
        void bind();
        void unbind();
    };

    class Camera {

    };
}

#endif /* GLL_HPP */