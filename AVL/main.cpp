#include <fstream>
#include <iostream>
// main currently configured to run solely the new AVLTree.hpp
#include "AVLTree.hpp"

int main() {
  std::ifstream input("input.txt");
  if (!input) {
    std::cerr << "Unable to open input.txt\n";
    return 1;
  }

  AVLTree<int> tree;
  int key;
  char option;

  while (input >> option) {
    if (option == 'I' || option == 'i') {
      if (!(input >> key)) {
        std::cerr << "Expected a key after I/i in input.txt\n";
        return 1;
      }
      tree.insert(key);
    } else if (option == 'D' || option == 'd') {
      if (!(input >> key)) {
        std::cerr << "Expected a key after D/d in input.txt\n";
        return 1;
      }
      tree.remove(key);
    } else {
      std::cerr << "Unknown command '" << option << "' in input.txt\n";
      return 1;
    }
  }

  if (!input.eof()) {
    std::cerr << "Error reading input.txt\n";
    return 1;
  }

  tree.display();
  return 0;
}