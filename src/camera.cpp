
#include "../include/camera.hpp"
#include "../include/gl.hpp"

void GL::Camera::drawTrackball(GL::ShapeRenderer& _shape){
    //_shape.drawCircle(Vec3(0.f), Vec3(1.f, 0.f, 0.f), 490.f, 500.f, RED);
    //_shape.drawCircle(Vec3(0.f), Vec3(0.f, 1.f, 0.f), 490.f, 500.f, GREEN);
    //_shape.drawCircle(Vec3(0.f), Vec3(0.f, 0.f, 1.f), 490.f, 500.f, BLUE);
    //_shape.drawSphere(Vec3(-100.f), 300.f, RED);
    //_shape.drawSphere(Vec3(0.f), 150.f, BLUE);
    //_shape.drawSphere(Vec3(100.f), 300.f, GREEN);

    //_shape.drawZylinder(Vec3(-250, 0.f, 0.f), Vec3(250.f, 0.f, 0.f), 50.f, RED);
    //_shape.drawLine(Vec3(0.f, -250.f, 0.f), Vec3(0.f, 250.f, 0.f), 50.f, GREEN);

    const Mat3 rotR = Mat3(glm::rotate(glm::radians(90.f), Vec3(0.f, 0.f, 1.f)));
    const Mat3 rotG = Mat3(glm::rotate(glm::radians(90.f), Vec3(0.f, 1.f, 0.f)));
    const Mat3 rotB = Mat3(glm::rotate(glm::radians(90.f), Vec3(1.f, 0.f, 0.f)));

    _shape.drawWheel(Vec3(0.f), 3.f, 600.f, 10.f, RED, rotR);
    _shape.drawWheel(Vec3(0.f), 3.f, 595.f, 10.f, GREEN, rotG);
    _shape.drawWheel(Vec3(0.f), 3.f, 590.f, 10.f, BLUE, rotB);
    
}

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