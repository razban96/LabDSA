cmake_minimum_required(VERSION 3.28)
project(LabDSA CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic -Wshadow -Wconversion)
endif()

# Target 1: AVL Tree Executable
add_executable(avl_tree
        AVL/main.cpp
        AVL/AVLTree.hpp
)

# Target 2: Heaps Executable
add_executable(heaps
        Heaps/heaps.cpp
)