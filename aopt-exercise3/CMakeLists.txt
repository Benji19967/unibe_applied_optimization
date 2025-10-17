cmake_minimum_required(VERSION 3.11)

PROJECT(AOPT_EXERCISES)

include(FetchContent)

FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)

if (WIN32)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build")

    # avoid linking errors, cf https://stackoverflow.com/questions/12540970/how-to-make-gtest-build-mdd-instead-of-mtd-by-default-using-cmake
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
else ()
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/bin")
endif ()


FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
  FetchContent_Populate(googletest)
  add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
endif()

enable_testing()
include(GoogleTest)


set(AOPT_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")


if(NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif()

# add library
add_library (AOPT INTERFACE)

add_library (AOPT::AOPT ALIAS AOPT)

#Set target properties
set(EIGEN3_INCLUDE_DIR $<BUILD_INTERFACE:${AOPT_EXTERNAL}/eigen>)

target_include_directories(AOPT
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        ${EIGEN3_INCLUDE_DIR}
)

add_subdirectory(EigenTutorial)
add_subdirectory(GridSearch)
add_subdirectory(CsvExporter)
add_subdirectory(ConvexityTests)
add_subdirectory(MassSpringSystem)
add_subdirectory(MassSpringProblemEvaluation)

