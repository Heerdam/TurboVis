#ifndef DEFINES_HPP
#define DEFINES_HPP

#define GLFW_INCLUDE_NONE
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
#include <utility>
#include <optional>
#include <array>
#include <random>
#include <algorithm>
#include <iostream>
#include <complex>
#include <future>

#define _USE_MATH_DEFINES
#include <math.h>

#define _SILENCE_CXX17_ADAPTOR_TYPEDEFS_DEPRECATION_WARNING
#define _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS
#include <Eigen/Eigen>
#include <Eigen/Core>

#include <immintrin.h>
#include <spdlog/spdlog.h>
#include <imgui.h>

#define GLM_ENABLE_EXPERIMENTAL
//#define GLM_FORCE_DEFAULT_ALIGNED_GENTYPES
#include <glm/ext/matrix_float4x4.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/ext.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/transform.hpp>

#include "allocator.hpp"
#include "exception.hpp"

using Vec2i = glm::ivec2;
using Vec3i = glm::ivec3;
using Vec4i = glm::ivec4;

using Vec2 = glm::vec2;
using Vec3 = glm::vec3;
using Vec4 = glm::vec4;

using Mat3 = glm::mat3x3;
using Mat4 = glm::mat4x4;

using Quat = glm::quat;

// ------------------- COLORS -------------------
const Vec4 RED = Vec4(1.f, 0.f, 0.f, 1.f);
const Vec4 GREEN = Vec4(0.f, 1.f, 0.f, 1.f);
const Vec4 BLUE = Vec4(0.f, 0.f, 1.f, 1.f);


#endif /* DEFINES_HPP */
