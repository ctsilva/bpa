cmake_minimum_required(VERSION 3.16)

project(bpa CXX)

# glm is header-only: use its CMake package when there is one, else the headers alone
find_package(glm CONFIG QUIET)
if(NOT TARGET glm::glm)
	find_path(GLM_INCLUDE_DIR glm/glm.hpp REQUIRED)
	add_library(glm::glm INTERFACE IMPORTED)
	target_include_directories(glm::glm INTERFACE ${GLM_INCLUDE_DIR})
endif()

# library
add_library(libbpa src/lib/bpa.cpp src/lib/bpa.h src/lib/IO.h)
target_link_libraries(libbpa PUBLIC glm::glm)
target_include_directories(libbpa PUBLIC src/lib)
target_compile_features(libbpa PUBLIC cxx_std_20)
target_compile_definitions(libbpa PUBLIC GLM_ENABLE_EXPERIMENTAL)
if(MSVC)
	target_compile_options(libbpa PRIVATE /W4)
else()
	target_compile_options(libbpa PRIVATE -Wall -Wextra -Wno-missing-field-initializers)
endif()
set_target_properties(libbpa PROPERTIES PREFIX "")

# driver
add_executable(bpa src/driver/main.cpp)
target_link_libraries(bpa PUBLIC libbpa)

# tests
include(CTest)
if(BUILD_TESTING)
	find_package(Catch2 3.0.1 QUIET)
endif()
if(NOT TARGET Catch2::Catch2WithMain)
	message(STATUS "Catch2 v3 not found: the tests are not built")
else()
	add_executable(tests test/computeBallCenter.cpp test/reconstruct.cpp)
	target_compile_definitions(tests PRIVATE TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/test/data" TEST_OUTPUT_DIR="${CMAKE_CURRENT_BINARY_DIR}")
	target_link_libraries(tests PRIVATE libbpa Catch2::Catch2WithMain)
	include(Catch)
	catch_discover_tests(tests)
endif()
