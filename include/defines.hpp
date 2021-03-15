#ifndef DEFINES_HPP
#define DEFINES_HPP

#include <vector>
#include <memory>
#include <cassert>
#include <cstring>

#include <Eigen/Eigen>
#include <immintrin.h>

using GL_Matrix4 = Eigen::Matrix<float, 4, 4, Eigen::ColMajor>;
using GL_Frustum = Eigen::Matrix<float, 6, 4>;

#endif /* DEFINES_HPP */