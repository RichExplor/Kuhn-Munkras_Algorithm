cmake_minimum_required(VERSION 2.8)

project(KM)

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -g")

find_package(Eigen3 REQUIRED)
find_package(OpenCV REQUIRED)
find_package(Ceres REQUIRED)

include_directories(
    ${OpenCV_INCLUDE_DIRS}
    ${CERES_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
)

add_executable(main main.cpp kuhn_munkras_algorithm.cpp)
target_link_libraries(main
    ${OpenCV_LIBS}
    ${EIGEN3_LIBS}
    ${CERES_LIBS}
    -pthread)