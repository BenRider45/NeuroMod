cmake_minimum_required(VERSION 3.10...3.50)
project(BCNNeuroModel VERSION 0.1)


set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include(FetchContent)

option(BUILD_TESTS "Build unit tests" ON)


# Gathering source files
file(GLOB_RECURSE LIB_SOURCES "src/BCNNeuroModel/*.cpp")
file(GLOB_RECURSE LIB_HEADERS "src/BCNNeuroModel/*.h" "src/BCNNeuroModel/*.hpp")

list(FILTER LIB_SOURCES EXCLUDE REGEX ".*main\\.cpp$")

add_library(BCNNeuroModel ${LIB_SOURCES} ${LIB_HEADERS})

target_compile_options(BCNNeuroModel PUBLIC -Wno-deprecated-copy)

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/BCNNeuroModel/main.cpp")
  add_executable(BCNNeuroModelDemo src/BCNNeuroModel/main.cpp)
    target_link_libraries(BCNNeuroModelDemo BCNNeuroModel)
    set_target_properties(BCNNeuroModelDemo PROPERTIES OUTPUT_NAME neuromod)
endif()

target_include_directories(BCNNeuroModel
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src/BCNNeuroModel
)


# Eigen dependency - check for local copy first, then find_package, then fetch
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/src/Eigen")
    message(STATUS "Using local Eigen copy from src/Eigen")
    target_include_directories(BCNNeuroModel PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src/Eigen")
else()
    find_package(Eigen3 QUIET)
    if(Eigen3_FOUND)
        message(STATUS "Found system Eigen3")
        target_link_libraries(BCNNeuroModel PUBLIC Eigen3::Eigen)
    else()
        message(STATUS "Eigen3 not found locally or on system, fetching from GitLab...")
        FetchContent_Declare(
            eigen
            GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
            GIT_TAG 3.4.0
            GIT_SHALLOW TRUE
        )
        FetchContent_MakeAvailable(eigen)
        target_link_libraries(BCNNeuroModel PUBLIC Eigen3::Eigen)
    endif()
endif()

target_compile_options(BCNNeuroModel PRIVATE
    $<$<CXX_COMPILER_ID:GNU>:-Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:Clang>:-Wall -Wextra -Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:/W4>
)

# gflags dependency
find_package(gflags REQUIRED)
target_link_libraries(BCNNeuroModel PUBLIC gflags)

set_target_properties(BCNNeuroModel PROPERTIES POSITION_INDEPENDENT_CODE ON)


if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

install(TARGETS BCNNeuroModel
    EXPORT BCNNeuroModelTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
)

install(FILES ${LIB_HEADERS} DESTINATION include/BCNNeuroModel)

install(EXPORT BCNNeuroModelTargets
    FILE BCNNeuroModelTargets.cmake
    NAMESPACE BCNNeuroModel::
    DESTINATION lib/cmake/BCNNeuroModel
)

# Create config file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "BCNNeuroModelConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/BCNNeuroModelConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/BCNNeuroModelConfig.cmake"
    INSTALL_DESTINATION lib/cmake/BCNNeuroModel
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/BCNNeuroModelConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/BCNNeuroModelConfigVersion.cmake"
    DESTINATION lib/cmake/BCNNeuroModel
)
