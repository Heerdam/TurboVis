
#include "../include/camera.hpp"
#include "../include/gl.hpp"

/*
    code taken and adapted from https://github.com/Pascal-So/turbotrack
    with permission by Pascal Sommer, 2021
*/
Vec3 GL::Camera::shoemake_projection(const Vec2& _mousePos, float _radius) {
    const float r2 = _radius * _radius;
    const float d2 = glm::dot(_mousePos, _mousePos);

    if (d2 <= r2) {
        // sphere
        return {_mousePos[0], _mousePos[1], std::sqrt(r2 - d2)};
    } else {
        // scaled sphere
        const float factor = _radius / std::sqrt(d2);
        return {factor * _mousePos[0], factor *  _mousePos[1], 0};
    }
}

/*
    code taken and adapted from https://github.com/Pascal-So/turbotrack
    with permission by Pascal Sommer, 2021
*/
Vec3 GL::Camera::holroyd_projection(const Vec2& _mousePos, float _radius) {
    const float r2 = _radius * _radius;
    const float d2 = glm::dot(_mousePos, _mousePos);

    if (d2 <= r2 / 2) {
        // sphere
        return {_mousePos[0], _mousePos[1], std::sqrt(r2 - d2)};
    } else {
        // hyperbola
        return {_mousePos[0], _mousePos[1], r2 / 2 / std::sqrt(d2)};
    }
}

/*
    code taken and adapted from https://github.com/Pascal-So/turbotrack
    with permission by Pascal Sommer, 2021
*/
Quat GL::Camera::trackball_shoemake(const Vec2& _oldPos, const Vec2& _newPos, float _radius){
    const Vec3 p1 = glm::normalize(shoemake_projection(_oldPos, _radius));
    const Vec3 p2 = glm::normalize(shoemake_projection(_newPos, _radius));
    return glm::rotation(p1, p2);
}

/*
    code taken and adapted from https://github.com/Pascal-So/turbotrack
    with permission by Pascal Sommer, 2021
*/
Quat GL::Camera::trackball_holroyd(const Vec2& _oldPos, const Vec2& _newPos, float _radius){
    const Vec3 p1 = glm::normalize(holroyd_projection(_oldPos, _radius));
    const Vec3 p2 = glm::normalize(holroyd_projection(_newPos, _radius));
    return glm::rotation(p1, p2);
}