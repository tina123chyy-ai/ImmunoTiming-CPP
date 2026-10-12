cmake_minimum_required(VERSION 3.16)
project(ImmunoTimingCPP VERSION 1.0.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()
add_library(immunotiming
    src/TumorModel.cpp
    src/TreatmentScheduler.cpp
    src/ODESolver.cpp
    src/ScheduleRanker.cpp
    src/Output.cpp)
target_include_directories(immunotiming PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
add_executable(immunotiming_cli src/main.cpp)
target_link_libraries(immunotiming_cli PRIVATE immunotiming)
enable_testing()
add_executable(immunotiming_tests tests/test_model.cpp)
target_link_libraries(immunotiming_tests PRIVATE immunotiming)
add_test(NAME immunotiming_tests COMMAND immunotiming_tests)
