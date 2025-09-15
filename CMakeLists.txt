cmake_minimum_required(VERSION 3.10)
project(breast_cancer)

set(CMAKE_CXX_STANDARD 17)

## This is the Apple compilation line
## set(CMAKE_CXX_FLAGS "-std=c++14 -g -O0 ${OpenMP_CXX_FLAGS}")

## This is the Windows compilation line
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "-fopenmp -std=c++17 -O3")

## Enable AddressSanitizer -- disabling for now
## set(ASAN_FLAGS "-fsanitize=address -fno-omit-frame-pointer -g")
## set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ASAN_FLAGS}")
## set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} ${ASAN_FLAGS}")

include_directories(inc)
file(GLOB SOURCES "src/*.cpp")
add_executable(breast_cancer ${SOURCES})

target_compile_options(breast_cancer PRIVATE -g -O0)
## Adding compile options for Address Sanitizer in Windows -- removing for now
## target_compile_options(breast_cancer PRIVATE -fsanitize=address)
## target_link_directories(breast_cancer PRIVATE "$ENV{ProgramFiles\(x86\)}/Microsoft Visual Studio/2019/Professional/VC/Tools/Llvm/x64/lib/clang/10.0.0/lib/windows")
## target_link_libraries(breast_cancer PRIVATE clang_rt.asan_dynamic-x86_64 clang_rt.asan_dynamic_runtime_thunk-x86_64)
## target_link_options(breast_cancer PRIVATE /wholearchive:clang_rt.asan_dynamic_runtime_thunk-x86_64.lib)


## OpenMP calls
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(breast_cancer PUBLIC OpenMP::OpenMP_CXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${OpenMP_SHARED_LINKER_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

## These are for Apple compilation
## set(CMAKE_CXX_FLAGS "-std=c++14 -g -O1 -fsanitize=address -fno-omit-frame-pointer")
## set(OpenMP_CXX_FLAGS "-Xclang -fopenmp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp")
## set(OpenMP_CXX_LIB_NAMES "omp")
## set(OpenMP_omp_LIBRARY omp)