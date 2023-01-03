
add_library( turbodorn INTERFACE )
target_include_directories(turbodorn INTERFACE "${CMAKE_CURRENT_LIST_DIR}/include")
set(LIBS_TT
  "nlohmann_json::nlohmann_json"
  "tsl::robin_map"
  "libmorton"
  "HDF5::HDF5"
  "OpenMP::OpenMP_CXX"
  "HighFive"
  "lodepng"
)
target_include_directories(turbodorn INTERFACE "${eigen3_SOURCE_DIR}")
target_link_libraries(turbodorn INTERFACE "${LIBS_TT}")

# ----------- Tests -----------
add_executable(turbodorn_test "${CMAKE_CURRENT_LIST_DIR}/tests/test_turbodorn.cpp") 
target_link_libraries(turbodorn_test turbodorn Catch2)