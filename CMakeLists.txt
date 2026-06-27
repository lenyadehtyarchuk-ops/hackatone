cmake_minimum_required(VERSION 3.16)
project(trn CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native")

find_package(OpenCV 4 REQUIRED COMPONENTS core imgproc imgcodecs highgui video)
find_package(GDAL REQUIRED)

find_package(OpenMP)

add_executable(trn
    src/main.cpp
    src/nmea_parser.cpp
    src/dem_loader.cpp
    src/terrain_correlator.cpp
    src/kalman_filter.cpp
    src/particle_filter.cpp
    src/keypoint_db.cpp
    src/contour_pf.cpp
    src/ridge_map.cpp
    src/line_navigator.cpp
    src/visualizer.cpp
)

target_include_directories(trn PRIVATE
    ${OpenCV_INCLUDE_DIRS}
    ${GDAL_INCLUDE_DIRS}
    src/
)

target_link_libraries(trn PRIVATE
    ${OpenCV_LIBS}
    ${GDAL_LIBRARIES}
)

if(OpenMP_CXX_FOUND)
    target_link_libraries(trn PRIVATE OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP enabled")
endif()

target_compile_options(trn PRIVATE -Wall -Wextra -Wpedantic)
