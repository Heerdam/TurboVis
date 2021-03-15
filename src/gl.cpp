#include "../include/gl.hpp"

GL_Frustum GL::Frustum::extractPlanes(const GL_Matrix4& _m){  

    GL_Frustum out;

    const __m128 _1x = _mm_load_ps(&_m(0, 0));
    const __m128 _2x = _mm_load_ps(&_m(1, 0));
    const __m128 _3x = _mm_load_ps(&_m(2, 0));
    const __m128 _4x = _mm_load_ps(&_m(3, 0));
    
    //plane 0
    {
        const __m128 p = _mm_add_ps(_4x, _1x);
        std::memcpy(out.data(), &p, sizeof(float) * 4);
    }

    //plane 1
    {
        const __m128 p = _mm_sub_ps(_4x, _1x);
        std::memcpy(out.data() + 4, &p, sizeof(float) * 4);
    }

    //plane 2
    {
        const __m128 p = _mm_sub_ps(_4x, _2x);
        std::memcpy(out.data() + 2*4, &p, sizeof(float) * 4);
    }

    //plane 3
    {
        const __m128 p = _mm_add_ps(_4x, _2x);
        std::memcpy(out.data() + 3*4, &p, sizeof(float) * 4);
    }

    //plane 4
    {
        const __m128 p = _mm_add_ps(_4x, _3x);
        std::memcpy(out.data() + 4*4, &p, sizeof(float) * 4);
    }

    //plane 5
    {
        const __m128 p = _mm_sub_ps(_4x, _3x);
        std::memcpy(out.data() + 5*4, &p, sizeof(float) * 4);
    }

    out.rowwise().normalize();
    return out;
}

bool GL::Frustum::isPointInside(const GL_Frustum& _planes, const Eigen::Vector3f& _point){
    for(size_t i = 0; i < 6; ++i){
        const float d = _planes.row(i).head(3).dot(_point) + _planes(i, 3);
        if(d <= 0.f)
            return false;
    }
    return true;
}