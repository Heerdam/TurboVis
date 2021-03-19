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
#include <filesystem>
#include <chrono>

#include <Eigen/Eigen>

#include <immintrin.h>
#include <spdlog/spdlog.h>
#include <imgui.h>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/ext/matrix_float4x4.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/ext.hpp>

using GL_Matrix4 = Eigen::Matrix<float, 4, 4, Eigen::ColMajor>;
using GL_Frustum = Eigen::Matrix<float, 6, 4>;


#endif /* DEFINES_HPP */