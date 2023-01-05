cmake_minimum_required(VERSION 3.16)

project(TurboVis)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g" )
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native")

include(FetchContent)

add_executable(tv "src/main.cpp")
add_executable(cam "src/camera_tool.cpp")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

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
	#target_include_directories(tv PUBLIC "${eigen3_SOURCE_DIR}")
	#target_include_directories(cam PUBLIC "${eigen3_SOURCE_DIR}")
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

# ------------------ Json ------------------
FetchContent_Declare(json
    GIT_REPOSITORY      https://github.com/ArthurSonzogni/nlohmann_json_cmake_fetchcontent.git
	GIT_TAG				v3.11.2
    GIT_SHALLOW         TRUE
)
FetchContent_MakeAvailable(json)

# ------------------ ROBIN-MAP ------------------
message(STATUS "Fetching robin map...")
FetchContent_Declare(robin_map
  GIT_REPOSITORY        https://github.com/Tessil/robin-map.git
  GIT_TAG               v1.0.1
  GIT_SHALLOW           TRUE
)
FetchContent_MakeAvailable(robin_map)

# ------------------ LIBMORTON ------------------
message(STATUS "Fetching libmorton...")
FetchContent_Declare(
    libmorton 
    GIT_REPOSITORY      https://github.com/Forceflow/libmorton.git
    GIT_TAG             v0.2.10
    GIT_SHALLOW         TRUE
)
FetchContent_MakeAvailable(libmorton)

# ------------------ CATCH2 ------------------
message(STATUS "Fetching Catch2...")
FetchContent_Declare(
		Catch2 
		GIT_REPOSITORY      https://github.com/catchorg/Catch2.git
		GIT_TAG             v2.13.6
		GIT_SHALLOW	        TRUE
)

find_package(HDF5 REQUIRED)
find_package(OpenMP REQUIRED)

option(CATCH_INSTALL_DOCS "Install documentation alongside library" OFF)
option(CATCH_INSTALL_EXTRAS "Install extras alongside library" OFF)
FetchContent_MakeAvailable(Catch2)

include("${CMAKE_CURRENT_SOURCE_DIR}/module/TurboDorn/module.cmake")

# ------------------------------------
target_link_libraries(tv turbodorn)
target_link_libraries(cam turbodorn)
