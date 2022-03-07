cmake_minimum_required(VERSION 3.16)

project(TurboVis)

if (MSVC)
	set(CMAKE_GENERATOR_PLATFORM x64)
endif()	
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if (MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj /openmp")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /Od ")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2 /DNDEBUG")
else()
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-rtti -fno-exceptions -mavx2")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG -O0 -g3 -ggdb -fsanitize=address -fsanitize=undefined -D_GLIBCXX_DEBUG" )
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -ffast-math -DNDEBUG -flto -fomit-frame-pointer")
endif()

include(FetchContent)

add_executable(tv "src/main.cpp")

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
	target_include_directories(tv PUBLIC "${eigen3_SOURCE_DIR}")
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

# ------------------ JSON ------------------
message(STATUS "Fetching json...")
FetchContent_Declare(json
  GIT_REPOSITORY https://github.com/ArthurSonzogni/nlohmann_json_cmake_fetchcontent.git
  GIT_TAG v3.10.4
  GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(json)

# ------------------ ROBIN-MAP ------------------
message(STATUS "Fetching robinmap...")
FetchContent_Declare(robin_map
  GIT_REPOSITORY https://github.com/Tessil/robin-map.git
  GIT_TAG v0.6.3
  GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(robin_map)

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

include("${CMAKE_CURRENT_SOURCE_DIR}/module/TurboDorn/module.cmake")

# ------------------------------------
target_link_libraries(tv turbodorn lodepng json robin_map)
