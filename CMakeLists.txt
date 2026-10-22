cmake_minimum_required(VERSION 3.28)
project(LabDSA CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic -Wshadow -Wconversion)
endif ()

# Target 1: AVL Tree Executable
add_executable(avl_tree
        AVL/main.cppa
        AVL/AVLTree.hpp
)
target_include_directories(avl_tree PRIVATE AVL)

# Target 2: Heaps Executable
add_executable(heaps
        Heaps/heaps.cpp
)
target_include_directories(heaps PRIVATE Heaps)

# Target 3: Trie Trees
add_executable(trie
        TrieTrees/trie.cpp)
target_include_directories(trie PRIVATE TrieTrees)