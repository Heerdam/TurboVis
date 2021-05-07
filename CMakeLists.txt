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
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g3 -ggdb -fsanitize=address -fsanitize=undefined -D_GLIBCXX_DEBUG" )
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Os -march=native -DNDEBUG -flto -fomit-frame-pointer")
	set(CMAKE_CXX_FLAGS_SCHRENZSPEED "-Ofast -funroll-loops -s -DNDEBUG -flto -fomit-frame-pointer -march=native")
endif()

include(FetchContent)

add_executable(tv
	"src/main.cpp"
	"src/data.cpp"
	"src/gl.cpp"
	"src/gui.cpp"
	"src/math.cpp"
	"src/util.cpp"
	"src/camera.cpp"
)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

# ------------------ Eigen ------------------
find_package (Eigen3 3.4.99 QUIET REQUIRED) 
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

	message(STATUS "Preparing lodepng...")

	set(LODEPNG_FILES
		"${lodepng_SOURCE_DIR}/lodepng.cpp"	
	)

	add_library(lodepng "${LODEPNG_FILES}")

	target_include_directories(tv PUBLIC "${lodepng_SOURCE_DIR}")
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

# ------------------ GLM ------------------
FetchContent_Declare(
        glm
        GIT_REPOSITORY https://github.com/g-truc/glm
		GIT_SHALLOW TRUE
		GIT_TAG 0.9.9.8
)
FetchContent_GetProperties(glm)
FetchContent_MakeAvailable(glm)

# ------------------ OMP ------------------
#find_package(OpenMP REQUIRED)

# ------------------ MRPT ------------------
message(STATUS "Fetching mrpt...")
FetchContent_Declare(
	mrpt
	GIT_REPOSITORY https://github.com/vioshyvo/mrpt.git
	GIT_SHALLOW TRUE
)
FetchContent_GetProperties(mrpt)

if(NOT mrpt_POPULATED)
	FetchContent_Populate(mrpt)

	message(STATUS "Preparing mrpt...")
	target_include_directories(tv PUBLIC "${mrpt_SOURCE_DIR}/cpp")
endif()

# ------------------------------------

set(TV_LIBS 
		"glfw"	
		"spdlog::spdlog"	
		"glad"
		"imgui"
		"lodepng"
		"libmorton" 
		"glm"
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

# ------------------ GnuPlot ------------------
	message(STATUS "Fetching gnuplot...")
	FetchContent_Declare(
		gp
		GIT_REPOSITORY https://github.com/sciplot/sciplot.git
		GIT_SHALLOW TRUE
		GIT_TAG v0.2.2
	)

	option(SCIPLOT_BUILD_EXAMPLES "Build examples" OFF)
	option(SCIPLOT_BUILD_TESTS "Build tests" OFF)
	option(SCIPLOT_BUILD_DOCS "Build documentation" OFF)

	FetchContent_GetProperties(gp)
	if(NOT gp_POPULATED)
		FetchContent_Populate(gp)
		add_subdirectory("${gp_SOURCE_DIR}" "${gp_BINARY_DIR}")
	endif()	

	# ------------------ Targets ------------------
	add_executable(VQ tests/VQ/main.cpp	)
	add_executable(SHAPES tests/SHAPES/main.cpp	)
	add_executable(PCA_STATIC tests/PCA_STATIC/main.cpp	)
	add_executable(VPM tests/VPM/main.cpp	)

	#VQ
	target_include_directories(VQ PUBLIC "${EIGEN3_INCLUDE_DIR}")
	target_include_directories(VQ PUBLIC "${imgui_SOURCE_DIR}")

	#SHAPES
	target_include_directories(SHAPES PUBLIC "${EIGEN3_INCLUDE_DIR}")
	target_include_directories(SHAPES PUBLIC "${imgui_SOURCE_DIR}")
	target_include_directories(SHAPES PUBLIC "${lodepng_SOURCE_DIR}")

	#PCA_STATIC
	target_include_directories(PCA_STATIC PUBLIC "${EIGEN3_INCLUDE_DIR}")
	target_include_directories(PCA_STATIC PUBLIC "${imgui_SOURCE_DIR}")

	#VPM
	target_include_directories(VPM PUBLIC "${EIGEN3_INCLUDE_DIR}")
	target_include_directories(VPM PUBLIC "${imgui_SOURCE_DIR}")
	target_include_directories(VPM PUBLIC "${lodepng_SOURCE_DIR}")

	if (MSVC)
		target_link_libraries(VQ "${TV_LIBS}" sciplot )
		target_link_libraries(SHAPES "${TV_LIBS}" sciplot )
		target_link_libraries(PCA_STATIC "${TV_LIBS}" sciplot )
		target_link_libraries(VPM "${TV_LIBS}" sciplot )
	else()
		target_link_libraries(VQ "${TV_LIBS_LINUX}" sciplot )
		target_link_libraries(SHAPES "${TV_LIBS_LINUX}" sciplot )
		target_link_libraries(PCA_STATIC "${TV_LIBS_LINUX}" sciplot )
		target_link_libraries(VPM "${TV_LIBS_LINUX}" sciplot )
	endif()
#endif()