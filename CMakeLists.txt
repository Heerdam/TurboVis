cmake_minimum_required(VERSION 3.16)

project(TurboVis)

set(CMAKE_GENERATOR_PLATFORM x64)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if (MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj /openmp")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Od ")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2 /DNDEBUG")
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-rtti -mavx -mavx2 -fopenmp")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG -O0 -g3 -ggdb -fsanitize=address -fsanitize=undefined -D_GLIBCXX_DEBUG" )
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Os -march=native -DNDEBUG -flto -fomit-frame-pointer")
	set(CMAKE_CXX_FLAGS_SCHRENZSPEED "-Ofast -funroll-loops -s -DNDEBUG -flto -fomit-frame-pointer -march=native")
endif()

include(FetchContent)

set(TV_FILES	
	"src/gui.cpp"
	"src/util.cpp"
	"src/camera.cpp"
)

add_executable(tv
	"${TV_FILES}"
	"src/main.cpp"
)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

find_package (HDF5 REQUIRED)

# ------------------ Eigen ------------------
message(STATUS "Fetching Eigen3...")
FetchContent_Declare(
	eigen3
	GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
	GIT_SHALLOW TRUE
	GIT_TAG        3.4
)
FetchContent_GetProperties(eigen3)

if(NOT eigen3_POPULATED)
	FetchContent_Populate(eigen3)
	target_include_directories(tv PUBLIC "${eigen3_SOURCE_DIR}")
endif()

# ------------------ HIGHFIVE ------------------
message(STATUS "Fetching HighFive...")
set(USE_BOOST OFF CACHE BOOL "Enable Boost Support")
set(USE_EIGEN OFF CACHE BOOL "Enable Eigen testing")
set(USE_XTENSOR OFF CACHE BOOL "Enable xtensor testing")
set(USE_OPENCV OFF CACHE BOOL "Enable OpenCV testing")
FetchContent_Declare(
	HighFive 
	GIT_REPOSITORY https://github.com/BlueBrain/HighFive.git
	GIT_SHALLOW TRUE
	GIT_TAG v2.3
)
FetchContent_GetProperties(HighFive)
FetchContent_MakeAvailable(HighFive)
#if(NOT HighFive_POPULATED)
	#FetchContent_Populate(HighFive)
	#option(HIGHFIVE_UNIT_TESTS "Enable unit tests" OFF)
	#option(HIGHFIVE_EXAMPLES "Compile examples" OFF)
	##option(HIGHFIVE_PARALLEL_HDF5 "Enable Parallel HDF5 support" ON)
	#option(HIGHFIVE_BUILD_DOCS "Enable documentation building" OFF)
	#add_subdirectory("${HighFive_SOURCE_DIR}" "${HighFive_BINARY_DIR}")
#endif()

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

# ------------------ LODEPNG ------------------
message(STATUS "Fetching lodepng...")
FetchContent_Declare(
	lodepng
	GIT_REPOSITORY https://github.com/lvandeve/lodepng.git
	GIT_SHALLOW TRUE
)
FetchContent_GetProperties(lodepng)

if(NOT lodepng_POPULATED)
	FetchContent_Populate(lodepng)

	set(LODEPNG_FILES
		"${lodepng_SOURCE_DIR}/lodepng.cpp"	
	)

	add_library(lodepng "${LODEPNG_FILES}")

	target_include_directories(tv PUBLIC "${lodepng_SOURCE_DIR}")
endif()

# ------------------ GLM ------------------
message(STATUS "Fetching GLM...")
FetchContent_Declare(
        glm
        GIT_REPOSITORY https://github.com/g-truc/glm
		GIT_SHALLOW TRUE
		GIT_TAG 0.9.9.8
)
FetchContent_MakeAvailable(glm)

# ------------------------------------

set(TV_LIBS 
		"glfw"	
		"spdlog::spdlog"	
		"glad"
		"imgui"
		"lodepng"
		"glm"
		"HighFive"
		#"OpenMP::OpenMP_CXX"
)

set(TV_LIBS_LINUX "${TV_LIBS}")
list(APPEND TV_LIBS_LINUX "dl")

if (MSVC)
	target_link_libraries(tv "${TV_LIBS}" )
else()
	target_link_libraries(tv "${TV_LIBS_LINUX}"	)
endif()

#if(TURBOVIS_BUILD_TESTS)

	# ------------------ CATCH2 ------------------
	message(STATUS "Fetching Catch2...")
	FetchContent_Declare(
			Catch2 
			GIT_REPOSITORY      https://github.com/catchorg/Catch2.git
			GIT_TAG             v2.13.6
			GIT_SHALLOW	        TRUE
	)

	option(CATCH_INSTALL_DOCS "Install documentation alongside library" OFF)
	option(CATCH_INSTALL_EXTRAS "Install extras alongside library" OFF)
	FetchContent_MakeAvailable(Catch2)

	# ------------------ Target ------------------
	add_executable( tv_tests "tests/main.cpp" "${TV_FILES}" )

	target_include_directories( tv_tests PUBLIC "${eigen3_SOURCE_DIR}" )
	target_include_directories( tv_tests PUBLIC "${imgui_SOURCE_DIR}" )
	target_include_directories( tv_tests PUBLIC "${lodepng_SOURCE_DIR}" )

	if (MSVC)
		target_link_libraries( tv_tests "${TV_LIBS}" Catch2 )
	else()
		target_link_libraries( tv_tests "${TV_LIBS_LINUX}" Catch2 )
	endif()
#endif()
