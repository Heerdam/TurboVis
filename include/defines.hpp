#ifndef DEFINES_HPP
#define DEFINES_HPP

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <vector>
#include <memory>
#include <cassert>
#include <cstring>
#include <fstream>
#include <functional>
#include <future>
#include <unordered_map>
#include <tuple>

#include <Eigen/Eigen>
#include <immintrin.h>
#include <spdlog/spdlog.h>

using GL_Matrix4 = Eigen::Matrix<float, 4, 4, Eigen::ColMajor>;
using GL_Frustum = Eigen::Matrix<float, 6, 4>;


#endif /* DEFINES_HPP */