cmake_minimum_required(VERSION 3.11)

project(TurboVis)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if (MSVC)
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2 /DNDEBUG")
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-rtti -mavx -mavx2")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -ggdb" )
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -DNDEBUG -flto")
endif()

include(FetchContent)

add_executable(tv
	src/main.cpp
	src/data.cpp
	src/gl.cpp
	src/gui.cpp
)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

# ------------------ Eigen ------------------
find_package (Eigen3 3.3.9 QUIET REQUIRED) 
target_include_directories(tv PUBLIC "${EIGEN3_INCLUDE_DIR}")

# ------------------ GLFW ------------------
message(STATUS "Fetching GLFW...")
FetchContent_Declare(
	glfw
	GIT_REPOSITORY https://github.com/glfw/glfw 
	GIT_SHALLOW TRUE
)

FetchContent_GetProperties(glfw)
if(NOT glfw_POPULATED)
	FetchContent_Populate(glfw)

	set(GLFW_BUILD_EXAMPLES OFF CACHE INTERNAL "Build the GLFW example programs")
	set(GLFW_BUILD_TESTS OFF CACHE INTERNAL "Build the GLFW test programs")
	set(GLFW_BUILD_DOCS OFF CACHE INTERNAL "Build the GLFW documentation")
	set(GLFW_INSTALL OFF CACHE INTERNAL "Generate installation target")

	add_subdirectory(${glfw_SOURCE_DIR} ${glfw_BINARY_DIR})
endif()

# ------------------ SPDLOG ------------------
message(STATUS "Fetching SPDLOG...")
FetchContent_Declare(
	spdlog
	GIT_REPOSITORY https://github.com/gabime/spdlog.git
	GIT_TAG        v1.8.2
	GIT_SHALLOW TRUE
)
FetchContent_GetProperties(spdlog)
FetchContent_MakeAvailable(spdlog)

# ------------------ GLAD ------------------ 

set(GLAD_API "gl=4.3" CACHE STRING "API type/version pairs, like \"gl=4.3,gles=\", no version means latest")
set(GLAD_EXTENSIONS "GL_ARB_buffer_storage" CACHE STRING "Path to extensions file or comma separated list of extensions, if missing all extensions are included")
FetchContent_Declare(
	glad
	GIT_REPOSITORY https://github.com/Dav1dde/glad.git
	GIT_TAG        v0.1.34
	GIT_SHALLOW TRUE
)
FetchContent_GetProperties(glad)
if(NOT glad_POPULATED)
	FetchContent_Populate(glad)
	add_subdirectory("${glad_SOURCE_DIR}" "${glad_BINARY_DIR}")
endif()
 
# ------------------ IMGUI ------------------
message(STATUS "Fetching imgui...")
FetchContent_Declare(
	imgui
	GIT_REPOSITORY https://github.com/ocornut/imgui.git
	GIT_SHALLOW TRUE
)
FetchContent_GetProperties(imgui)

if(NOT imgui_POPULATED)
	FetchContent_Populate(imgui)

	message(STATUS "Preparing imgui...")

	set(IMGUI_FILES
		"${imgui_SOURCE_DIR}/imgui_demo.cpp"
		"${imgui_SOURCE_DIR}/imgui_draw.cpp"
		"${imgui_SOURCE_DIR}/imgui_widgets.cpp"
		"${imgui_SOURCE_DIR}/imgui_tables.cpp"
		"${imgui_SOURCE_DIR}/imgui.cpp"		
	)

	add_library(imgui "${IMGUI_FILES}")

	target_include_directories(tv PUBLIC "${imgui_SOURCE_DIR}")
endif()

# ------------------ LIBMORTON ------------------
message(STATUS "Fetching LIBMORTON...")
FetchContent_Declare(
	libmorton 
	GIT_REPOSITORY https://github.com/Forceflow/libmorton.git
	GIT_TAG        v0.2.7
	GIT_SHALLOW TRUE
)
option(BUILD_TESTING "Build unit tests for libmorton" OFF)
FetchContent_GetProperties(libmorton)
FetchContent_MakeAvailable(libmorton)


FetchContent_Declare(
        glm
        GIT_REPOSITORY https://github.com/g-truc/glm
		GIT_SHALLOW TRUE
		GIT_TAG 0.9.9.8
)
FetchContent_GetProperties(glm)
if(NOT glm_POPULATED)
	FetchContent_Populate(glm)
	add_subdirectory("${glm_SOURCE_DIR}" "${glm_BINARY_DIR}")
endif()

if (MSVC)
	target_link_libraries(tv
		glfw		
		spdlog::spdlog	
		glad
		imgui
		libmorton 
		glm
	)
else()
	target_link_libraries(tv
		glfw		
		spdlog::spdlog	
		glad
		imgui
		dl
		libmorton 
		glm
	)
endif()