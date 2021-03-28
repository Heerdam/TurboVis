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
#include <utility>
#include <optional>
#include <array>
#include <random>
#include <algorithm>

#include <Eigen/Eigen>
#include <Eigen/Core>

#include <immintrin.h>
#include <spdlog/spdlog.h>
#include <imgui.h>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/ext/matrix_float4x4.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/ext.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>

#include "allocator.hpp"

using Vec2i = glm::ivec2;
using Vec3i = glm::ivec3;
using Vec4i = glm::ivec4;

using Vec2 = glm::vec2;
using Vec3 = glm::vec3;
using Vec4 = glm::vec4;

using Mat3 = glm::mat3x3;
using Mat4 = glm::mat4x4;

using Quat = glm::quat;


#endif /* DEFINES_HPP */