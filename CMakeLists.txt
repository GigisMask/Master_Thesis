cmake_minimum_required(VERSION 3.0.0)
project(RandPot VERSION 0.1.0)

include(CTest)
enable_testing()

set(CMAKE_CXX_FLAGS "-lpthread -lgsl -lgslcblas -lfftw3 -lm")
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost 1.45.0 COMPONENTS filesystem) 

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    aux_source_directory(src src)
    add_executable(${PROJECT_NAME}.o main.cpp ${src})
    set_target_properties( ${PROJECT_NAME}.o PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})
    target_link_libraries(${PROJECT_NAME}.o ${Boost_LIBRARIES})
endif()



set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
