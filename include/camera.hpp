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

        static Vec3 shoemake_projection(const Vec2&, float);
        static Vec3 holroyd_projection(const Vec2&, float);
        static Quat trackball_shoemake(const Vec2&, const Vec2&, float);
        static Quat trackball_holroyd(const Vec2&, const Vec2&, float);

        [[nodiscard]] Vec3& operator*(Vec3& _in) const noexcept {
            _in = normalize(_in.x * right + _in.y * up + _in.z * dir);
            return _in;
        }
    };

}

#endif /* CAMERA_HPP */