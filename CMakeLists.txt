cmake_minimum_required(VERSION 3.0.0)
project(RandPot VERSION 0.1.0)

include(CTest)
enable_testing()

set(CMAKE_CXX_FLAGS "-lpthread -lgsl -lgslcblas -lfftw3 -lm")
aux_source_directory(src src)
add_executable(${PROJECT_NAME}.o main.cpp ${src})
set_target_properties( ${PROJECT_NAME}.o PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
