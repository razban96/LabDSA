#pragma once

#include <algorithm>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

template <typename T> class AVLTree {
private:
// contains variable data type, a height of min 1
// and unique pointers to left and right nodes
  struct Node {
    T data;
    int height{1};
    std::unique_ptr<Node> left{nullptr};
    std::unique_ptr<Node> right{nullptr};

    // blocks the node from taking any other form, transfers ownership of ptr with move
    explicit Node(T val) : data(std::move(val)) {}
  };
  // the root is initialized to nullptr
  std::unique_ptr<Node> root{nullptr};

  // Private helper declarations
  int getHeight(const Node *node) const;
  int getBalanceFactor(const Node *node) const;
  const Node *getMin(const Node *node) const;

  // rotations and rebalance
  std::unique_ptr<Node> rotateLeft(std::unique_ptr<Node> node);
  std::unique_ptr<Node> rotateRight(std::unique_ptr<Node> node);
  std::unique_ptr<Node> rebalance(std::unique_ptr<Node> node);

  // operations on nodes
  void insertNode(std::unique_ptr<Node> &node, T &&val);
  void deleteNode(std::unique_ptr<Node> &node, const T &val);
  std::unique_ptr<Node> extractMin(std::unique_ptr<Node> &node);
  void displayTree(const Node *node) const;

public:
    // default constructor
  AVLTree() = default;

  // Public API interface
  void insert(T val);
  void remove(const T &val);
  void display() const;
};

// ===================================================
// METHOD DEFINITIONS
// ===================================================

template <typename T> void AVLTree<T>::insert(T val) {
  insertNode(root, std::move(val));
}

template <typename T> void AVLTree<T>::remove(const T &val) {
  deleteNode(root, val);
}

template <typename T> void AVLTree<T>::display() const {
  displayTree(root.get());
}

template <typename T>
void AVLTree<T>::displayTree(const Node *node) const {
  if (node == nullptr) {
    return;
  }

  std::vector<std::vector<const Node *>> levels{{node}};
  while (true) {
    const auto &currentLevel = levels.back();
    std::vector<const Node *> nextLevel;
    nextLevel.reserve(currentLevel.size() * 2);
    bool hasChildren = false;

    for (const Node *current : currentLevel) {
      if (current == nullptr) {
        nextLevel.push_back(nullptr);
        nextLevel.push_back(nullptr);
        continue;
      }

      nextLevel.push_back(current->left.get());
      nextLevel.push_back(current->right.get());
      hasChildren = hasChildren || current->left != nullptr ||
                    current->right != nullptr;
    }

    if (!hasChildren) {
      break;
    }
    levels.push_back(std::move(nextLevel));
  }

  std::vector<std::vector<std::string>> labels;
  std::size_t labelWidth = 0;
  for (const auto &level : levels) {
    std::vector<std::string> levelLabels;
    levelLabels.reserve(level.size());
    for (const Node *current : level) {
      if (current == nullptr) {
        levelLabels.emplace_back();
        continue;
      }

      std::ostringstream label;
      label << current->data << " (bf=" << getBalanceFactor(current) << ")";
      levelLabels.push_back(label.str());
      labelWidth = std::max(labelWidth, levelLabels.back().size());
    }
    labels.push_back(std::move(levelLabels));
  }

  const std::size_t slotWidth = labelWidth + 2;
  const std::size_t canvasWidth = levels.back().size() * slotWidth;
  const auto printLine = [](const std::string &line) {
    const std::size_t last = line.find_last_not_of(' ');
    std::cout << (last == std::string::npos ? std::string() : line.substr(0, last + 1))
              << '\n';
  };

  for (std::size_t levelIndex = 0; levelIndex < levels.size(); ++levelIndex) {
    const auto &level = levels[levelIndex];
    const std::size_t slotSpan = canvasWidth / level.size();
    std::string nodeLine(canvasWidth, ' ');

    for (std::size_t index = 0; index < level.size(); ++index) {
      const std::string &label = labels[levelIndex][index];
      if (label.empty()) {
        continue;
      }

      const std::size_t center = index * slotSpan + slotSpan / 2;
      const std::size_t start = center - label.size() / 2;
      nodeLine.replace(start, label.size(), label);
    }
    printLine(nodeLine);

    if (levelIndex + 1 == levels.size()) {
      continue;
    }

    std::string branchLine(canvasWidth, ' ');
    for (std::size_t index = 0; index < level.size(); ++index) {
      const Node *current = level[index];
      if (current == nullptr) {
        continue;
      }

      const std::size_t center = index * slotSpan + slotSpan / 2;
      if (current->left != nullptr) {
        branchLine[center - slotSpan / 4] = '/';
      }
      if (current->right != nullptr) {
        branchLine[center + slotSpan / 4] = '\\';
      }
    }
    printLine(branchLine);
  }
}

template <typename T>
const typename AVLTree<T>::Node *AVLTree<T>::getMin(const Node *node) const {
  const Node *current = node;
  while (current != nullptr && current->left != nullptr) {
    current = current->left.get();
  }
  return current;
}

template <typename T> int AVLTree<T>::getHeight(const Node *node) const {
  return (node != nullptr) ? node->height : 0;
}

template <typename T> int AVLTree<T>::getBalanceFactor(const Node *node) const {
  if (node == nullptr) {
    return 0;
  }
  return getHeight(node->left.get()) - getHeight(node->right.get());
}

template <typename T>
std::unique_ptr<typename AVLTree<T>::Node>
AVLTree<T>::extractMin(std::unique_ptr<Node> &node) {
  if (node == nullptr) {
    return nullptr;
  }

  if (node->left == nullptr) {
    auto min = std::move(node);
    node = std::move(min->right);
    return min;
  }

  auto min = extractMin(node->left);
  node = rebalance(std::move(node));
  return min;
}

template <typename T>
std::unique_ptr<typename AVLTree<T>::Node>
AVLTree<T>::rotateRight(std::unique_ptr<Node> node) {
  if (node == nullptr || node->left == nullptr) {
    return node;
  }

  auto b = std::move(node->left);
  node->left = std::move(b->right);

  node->height =
      1 + std::max(getHeight(node->left.get()), getHeight(node->right.get()));

  b->right = std::move(node);

  b->height = 1 + std::max(getHeight(b->left.get()), getHeight(b->right.get()));
  return b;
}

template <typename T>
std::unique_ptr<typename AVLTree<T>::Node>
AVLTree<T>::rotateLeft(std::unique_ptr<Node> node) {
  if (node == nullptr || node->right == nullptr) {
    return node;
  }

  auto b = std::move(node->right);
  node->right = std::move(b->left);

  node->height =
      1 + std::max(getHeight(node->left.get()), getHeight(node->right.get()));

  b->left = std::move(node);

  b->height = 1 + std::max(getHeight(b->left.get()), getHeight(b->right.get()));
  return b;
}

template <typename T>
std::unique_ptr<typename AVLTree<T>::Node>
AVLTree<T>::rebalance(std::unique_ptr<Node> node) {
  if (node == nullptr)
    return nullptr;

  node->height =
      1 + std::max(getHeight(node->left.get()), getHeight(node->right.get()));

  int balance = getBalanceFactor(node.get());

  if (balance > 1) {
    if (getBalanceFactor(node->left.get()) < 0) {
      node->left = rotateLeft(std::move(node->left));
    }
    return rotateRight(std::move(node));
  }
  if (balance < -1) {
    if (getBalanceFactor(node->right.get()) > 0) {
      node->right = rotateRight(std::move(node->right));
    }
    return rotateLeft(std::move(node));
  }

  return node;
}

template <typename T>
void AVLTree<T>::insertNode(std::unique_ptr<Node> &node, T &&val) {
  if (node == nullptr) {
    node = std::make_unique<Node>(std::forward<T>(val));
    return;
  }

  if (val < node->data) {
    insertNode(node->left, std::forward<T>(val));
  } else if (node->data < val) {
    insertNode(node->right, std::forward<T>(val));
  } else {
    return;
  }

  node = rebalance(std::move(node));
}

template <typename T>
void AVLTree<T>::deleteNode(std::unique_ptr<Node> &node, const T &val) {
  if (node == nullptr) {
    return;
  }

  if (val < node->data) {
    deleteNode(node->left, val);
  } else if (node->data < val) {
    deleteNode(node->right, val);
  } else {
    // 0 or 1 child case
    if (node->left == nullptr) {
      node = std::move(node->right);
      return;
    }
    if (node->right == nullptr) {
      node = std::move(node->left);
      return;
    }

    // 2 children case
    auto successor = extractMin(node->right);
    successor->left = std::move(node->left);
    successor->right = std::move(node->right);
    node = rebalance(std::move(successor));
    return;
  }

  node = rebalance(std::move(node));
}