#ifndef GLL_HPP
#define GLL_HPP

#include "defines.hpp"

namespace GL {

    namespace Frustum {
        GL_Frustum extractPlanes(const GL_Matrix4&);
        bool isPointInside(const GL_Frustum&, const Eigen::Vector3f&);
    }

}

#endif /* GLL_HPP */