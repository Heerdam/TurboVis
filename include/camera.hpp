#ifndef CAMERA_HPP
#define CAMERA_HPP

#include "defines.hpp"

namespace GL {

    class Camera {
    public:
        const Vec3 upAxis = Vec3(0.f, 1.f, 0.f);

        Camera(int64_t _width, int64_t _height, float _fov, float _near, float _far) noexcept
            : width(_width), height(_height), fov(_fov), near(_near), far(_far) {}

        void update() noexcept {
            dir = normalize(target - position);
            right = normalize(cross(dir, upAxis));
            up = normalize(cross(right, dir));
        }

        bool hasMoved = true;

        Vec3 position, target;
        Vec3 up, right, dir;

        const float near, far, fov;
        const int64_t width, height;

        Mat4 combined;
        Mat4 comb;

        [[nodiscard]] static Vec3 shoemake_projection(const Vec2& /*_mousePos*/, float /*_radius*/) noexcept;
        [[nodiscard]] static Vec3 holroyd_projection(const Vec2& /*_mousePos*/, float /*_radius*/) noexcept;
        [[nodiscard]] static Quat trackball_shoemake(const Vec2& /*_oldPos*/, const Vec2& /*_newPos*/, float /*_radius*/) noexcept;
        [[nodiscard]] static Quat trackball_holroyd(const Vec2& /*_oldPos*/, const Vec2& /*_newPos*/, float /*_radius*/) noexcept;

        [[nodiscard]] Vec3& operator*(Vec3& _in) const noexcept {
            _in = glm::normalize(_in.x * right + _in.y * up + _in.z * dir);
            return _in;
        }
    };

}

inline Vec3 GL::Camera::shoemake_projection(const Vec2& _mousePos, float _radius) noexcept {
    const float r2 = _radius * _radius;
    const float d2 = glm::dot(_mousePos, _mousePos);

    if (d2 <= r2) {
        // sphere
        return { _mousePos[0], _mousePos[1], std::sqrt(r2 - d2) };
    } else {
        // scaled sphere
        const float factor = _radius / std::sqrt(d2);
        return { factor * _mousePos[0], factor *  _mousePos[1], 0} ;
    }
}

inline Vec3 GL::Camera::holroyd_projection(const Vec2& _mousePos, float _radius) noexcept {
    const float r2 = _radius * _radius;
    const float d2 = glm::dot(_mousePos, _mousePos);

    if (d2 <= r2 / 2) {
        // sphere
        return { _mousePos[0], _mousePos[1], std::sqrt(r2 - d2) };
    } else {
        // hyperbola
        return { _mousePos[0], _mousePos[1], r2 / 2 / std::sqrt(d2) };
    }
}

inline Quat GL::Camera::trackball_shoemake(const Vec2& _oldPos, const Vec2& _newPos, float _radius)noexcept {
    const Vec3 p1 = glm::normalize(shoemake_projection(_oldPos, _radius));
    const Vec3 p2 = glm::normalize(shoemake_projection(_newPos, _radius));
    return glm::rotation(p1, p2);
}

inline Quat GL::Camera::trackball_holroyd(const Vec2& _oldPos, const Vec2& _newPos, float _radius) noexcept{
    const Vec3 p1 = glm::normalize(holroyd_projection(_oldPos, _radius));
    const Vec3 p2 = glm::normalize(holroyd_projection(_newPos, _radius));
    return glm::rotation(p1, p2);
}

#endif /* CAMERA_HPP */