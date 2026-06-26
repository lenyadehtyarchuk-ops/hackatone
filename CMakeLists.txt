cmake_minimum_required(VERSION 3.18)
project(kazaki_car_racing LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release)
endif()

# ---------------- GDAL ----------------
find_package(GDAL QUIET)

if(NOT GDAL_FOUND)
    set(GDAL_DEV_ROOT "${CMAKE_SOURCE_DIR}/third_party/gdal-dev/usr")
    if(EXISTS "${GDAL_DEV_ROOT}/include/gdal/gdal_priv.h")
        set(GDAL_INCLUDE_DIR "${GDAL_DEV_ROOT}/include/gdal")
        set(_gdal_lib_candidates
            /usr/lib/x86_64-linux-gnu/libgdal.so
            /usr/lib/x86_64-linux-gnu/libgdal.so.38
        )
        set(GDAL_LIBRARY "")
        foreach(_candidate IN LISTS _gdal_lib_candidates)
            if(EXISTS "${_candidate}")
                set(GDAL_LIBRARY "${_candidate}")
                break()
            endif()
        endforeach()
        if(NOT GDAL_LIBRARY)
            message(FATAL_ERROR "libgdal not found in /usr/lib/x86_64-linux-gnu")
        endif()
        set(GDAL_LIBRARIES ${GDAL_LIBRARY})
        message(STATUS "Using vendored GDAL headers from ${GDAL_DEV_ROOT}")
    else()
        message(FATAL_ERROR
            "GDAL not found. Install dev package:\n"
            "  sudo apt-get install -y libgdal-dev\n"
            "Or run VS Code task: install deps (GDAL)")
    endif()
endif()

add_executable(gdal_test src/main.cpp)
target_include_directories(gdal_test PRIVATE ${GDAL_INCLUDE_DIR})
target_link_libraries(gdal_test PRIVATE ${GDAL_LIBRARIES})

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "GDAL include: ${GDAL_INCLUDE_DIR}")
message(STATUS "GDAL libraries: ${GDAL_LIBRARIES}")
