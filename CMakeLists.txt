cmake_minimum_required(VERSION 3.27)
project(PLR)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

add_library(PLR STATIC library.cpp)

enable_testing()

add_executable(
        PLRTest
        PLRTest.cc
)
target_link_libraries(
        PLRTest
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(PLRTest)
