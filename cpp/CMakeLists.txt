cmake_minimum_required(VERSION 3.10)

project(bezier)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-O3 -pthread")

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

# find_package(cdd REQUIRED)
# include_directories(/home/noel/repos/eigen-cddlib_install/include/eigen-cdd)

add_library(BEZIER SHARED
	src/Bezier.cpp
)

target_link_libraries(BEZIER Eigen3::Eigen)
# target_link_libraries(BEZIER cdd)

add_executable(main src/main.cpp)
target_link_libraries(main BEZIER)
add_library(cddeigen SHARED IMPORTED)
# set_target_properties(cddeigen PROPERTIES IMPORTED_LOCATION /home/noel/repos/eigen-cddlib_install/lib/libeigen-cdd.so)
# target_link_libraries(main cddeigen)

#################### GoogleTest ###########################
# include(FetchContent)
# FetchContent_Declare(
#   googletest
#   URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
# )
# FetchContent_MakeAvailable(googletest)
# enable_testing()
# add_executable(
#   hello_test
#   unittest/hello_test.cc
# )
# target_link_libraries(
#   hello_test
#   GTest::gtest_main
#   BEZIER
# )
# include(GoogleTest)
# gtest_discover_tests(hello_test)

############################### Pybind11 ##############################
# add_subdirectory(pybind11)
# pybind11_add_module(bezier "../bindings/bezier.cpp")
# target_link_libraries(bezier PUBLIC BEZIER)
# target_link_libraries(bezier PUBLIC cddeigen)
