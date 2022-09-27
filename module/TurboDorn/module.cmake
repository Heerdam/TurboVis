
add_library( turbodorn INTERFACE )
target_include_directories(turbodorn INTERFACE "${CMAKE_CURRENT_LIST_DIR}/include")
set(LIBS_TT
  #"nlohmann_json::nlohmann_json"
  #"tsl::robin_map"
)
target_include_directories(turbodorn INTERFACE "${eigen3_SOURCE_DIR}")
target_link_libraries(turbodorn INTERFACE)

# ----------- Tests -----------
add_executable(turbodorn_test "${CMAKE_CURRENT_LIST_DIR}/tests/test_turbodorn.cpp") 
target_link_libraries(turbodorn_test "${LIBS_TT}" Catch2 turbodorn)