cmake_minimum_required(VERSION 3.10)
project(breast_cancer)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "-std=c++14 -g -O0 ${OpenMP_CXX_FLAGS}")

## Enable AddressSanitizer
#set(ASAN_FLAGS "-fsanitize=address -fno-omit-frame-pointer -g")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ASAN_FLAGS}")
#set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} ${ASAN_FLAGS}")

include_directories(inc)
file(GLOB SOURCES "src/*.cpp")
add_executable(breast_cancer ${SOURCES})

target_compile_options(breast_cancer PRIVATE -g -O0)


find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(breast_cancer PUBLIC OpenMP::OpenMP_CXX)
endif()

#set(CMAKE_CXX_FLAGS "-std=c++14 -g -O1 -fsanitize=address -fno-omit-frame-pointer")
set(OpenMP_CXX_FLAGS "-Xclang -fopenmp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp")
set(OpenMP_CXX_LIB_NAMES "omp")
set(OpenMP_omp_LIBRARY omp)