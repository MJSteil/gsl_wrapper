cmake_minimum_required(VERSION 3.6)
project(gsl_wrapper)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE} -march=native -O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wno-unknown-pragmas -Wextra")

file(GLOB HEADERS *.hpp)
file(GLOB_RECURSE SOURCE source/*.cpp)

add_library(gsl_wrapper ${HEADERS} ${SOURCE})
target_link_libraries(gsl_wrapper gsl gslcblas)

set_target_properties(gsl_wrapper PROPERTIES
        SUFFIX -${CMAKE_BUILD_TYPE}.a
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)