// B-trees - DSA lab
//
// B-tree of order d (Bayer and McCreight, 1972), as in the examples of the course:
//  - every node holds at most 2d keys; an internal node with k keys has k + 1 children;
//  - every node except the root holds at least d keys;
//  - all leaves are on the same level;
//  - the keys of every node are in increasing order.
// With d = 2, a node holds 2 to 4 keys and has 3 to 5 children. In the definition of
// the course, where the order m is the maximum number of children, this is a B-tree
// of order m = 5.
//
// Build: g++ -std=c++11 -Wall -Wextra -o btree_menu_en btree_menu_en.cpp

#include <cerrno>
#include <climits>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

const int D = 2;          // the order of the tree
const int MAXK = 2 * D;   // the largest number of keys in a node
const int MINK = D;       // the smallest number of keys in a node other than the root

struct NodeB {
    int nr;                   // the number of keys in the node
    int key[MAXK + 2];        // key[1..nr]; key[MAXK + 1] briefly holds one extra key
    NodeB* child[MAXK + 2];   // child[0..nr]; child[i]: keys from key[i] to key[i+1]

    NodeB() : nr(0) {
        for (int i = 0; i < MAXK + 2; i++) { key[i] = 0; child[i] = nullptr; }
    }
    bool leaf() const { return child[0] == nullptr; }
};

class BTree {
public:
    BTree() : root(nullptr), count(0) {}
    ~BTree() { clear(); }
    BTree(const BTree&) = delete;
    BTree& operator=(const BTree&) = delete;

    int size() const { return count; }
    bool empty() const { return root == nullptr; }

    // searches for k; if path is not null, it receives the visited nodes, root first
    bool search(int k, std::vector<const NodeB*>* path = nullptr) const {
        const NodeB* n = root;
        while (n != nullptr) {
            if (path != nullptr) path->push_back(n);
            int i = position(n, k);
            if (i <= n->nr && n->key[i] == k) return true;
            n = n->child[i - 1];               // at a leaf, child[i - 1] is nullptr
        }
        return false;
    }

    // inserts k; returns false if the key is already there
    bool insert(int k) {
        if (root == nullptr) {
            root = new NodeB();
            root->nr = 1;
            root->key[1] = k;
            count = 1;
            return true;
        }
        int up = 0;
        NodeB* right = nullptr;
        Result r = insertIn(root, k, up, right);
        if (r == EXISTS) return false;
        if (r == SPLIT) {                      // the root was split: one more level
            NodeB* s = new NodeB();
            s->nr = 1;
            s->key[1] = up;
            s->child[0] = root;
            s->child[1] = right;
            root = s;
        }
        count++;
        return true;
    }

    // deletes k; returns false if the key is not there
    bool remove(int k) {
        if (root == nullptr || !removeFrom(root, k)) return false;
        count--;
        if (root->nr == 0) {                   // the root has no keys left
            NodeB* old = root;
            root = root->child[0];          // nullptr if the tree is now empty
            delete old;
        }
        return true;
    }

    // all the keys, in increasing order
    std::vector<int> keys() const {
        std::vector<int> out;
        inorder(root, out);
        return out;
    }

    // the keys in the interval [a, b], in increasing order
    std::vector<int> range(int a, int b) const {
        std::vector<int> out;
        inRange(root, a, b, out);
        return out;
    }

    int nodes() const { return countNodes(root); }

    // the height: the number of edges from the root to the leaves (-1 if empty)
    int height() const {
        int h = -1;
        for (const NodeB* n = root; n != nullptr; n = n->child[0]) h++;
        return h;
    }

    // checks the B-tree properties; why receives the first property broken
    bool check(std::string& why) const {
        if (root == nullptr) {
            if (count == 0) return true;
            why = "the tree is empty, but the key counter is not 0";
            return false;
        }
        int leafLevel = -1, n = 0;
        if (!checkNode(root, 0, true, false, 0, false, 0, leafLevel, n, why))
            return false;
        if (n != count) {
            why = "the number of keys does not match the counter";
            return false;
        }
        return true;
    }

    // prints the tree level by level; "|" separates children of different parents
    void print() const {
        if (root == nullptr) { std::cout << "  (the tree is empty)\n"; return; }
        typedef std::vector<const NodeB*> Group;   // the children of one parent
        std::vector<Group> level(1, Group(1, root));
        for (int d = 0; !level.empty(); d++) {
            std::vector<Group> next;
            std::cout << "  level " << d << ":";
            for (size_t g = 0; g < level.size(); g++) {
                if (g > 0) std::cout << " |";
                for (const NodeB* n : level[g]) {
                    std::cout << ' ' << str(n);
                    if (n->leaf()) continue;
                    next.push_back(Group());
                    for (int i = 0; i <= n->nr; i++) next.back().push_back(n->child[i]);
                }
            }
            std::cout << '\n';
            level.swap(next);
        }
    }

    void clear() {
        destroy(root);
        root = nullptr;
        count = 0;
    }

    // a node written as in the course, for example
    static std::string str(const NodeB* n) {
        std::string s = "[";
        for (int i = 1; i <= n->nr; i++) {
            if (i > 1) s += ' ';
            s += std::to_string(n->key[i]);
        }
        return s + "]";
    }

private:
    NodeB* root;
    int count;

    enum Result { EXISTS, DONE, SPLIT };

    // the first position i with k <= key[i], or nr + 1 if k exceeds every key
    static int position(const NodeB* n, int k) {
        int i = 1;
        while (i <= n->nr && k > n->key[i]) i++;
        return i;
    }

    // inserts k into the subtree rooted at n; if n is split, returns SPLIT,
    // with the median key in up and the new right node in right
    static Result insertIn(NodeB* n, int k, int& up, NodeB*& right) {
        int i = position(n, k);
        if (i <= n->nr && n->key[i] == k) return EXISTS;
        int x = k;                             // the key that goes into n
        NodeB* xr = nullptr;                   // its right child (null at leaves)
        if (!n->leaf()) {
            Result r = insertIn(n->child[i - 1], k, x, xr);
            if (r != SPLIT) return r;          // the child took the key without a split
        }
        for (int j = n->nr; j >= i; j--) {     // x goes to position i, xr to its right
            n->key[j + 1] = n->key[j];
            n->child[j + 1] = n->child[j];
        }
        n->key[i] = x;
        n->child[i] = xr;
        n->nr++;
        if (n->nr <= MAXK) return DONE;
        split(n, up, right);                   // n holds 2d + 1 keys: split it
        return SPLIT;
    }

    // n holds 2d + 1 keys: the first d stay in n, key d + 1 moves up to the parent,
    // and the last d move to the new node right
    static void split(NodeB* n, int& up, NodeB*& right) {
        right = new NodeB();
        right->child[0] = n->child[D + 1];
        n->child[D + 1] = nullptr;
        for (int j = 1; j <= D; j++) {
            right->key[j] = n->key[D + 1 + j];
            right->child[j] = n->child[D + 1 + j];
            n->child[D + 1 + j] = nullptr;
        }
        right->nr = D;
        up = n->key[D + 1];
        n->nr = D;
    }

    // deletes k from the subtree rooted at n; a child left with d - 1 keys is fixed
    static bool removeFrom(NodeB* n, int k) {
        int i = position(n, k);
        bool here = i <= n->nr && n->key[i] == k;
        if (n->leaf()) {
            if (!here) return false;
            for (int j = i; j < n->nr; j++) n->key[j] = n->key[j + 1];
            n->nr--;
            return true;
        }
        if (here) {
            // internal node: k is replaced with its predecessor, the largest key of
            // the left subtree; the predecessor is then deleted from its leaf
            const NodeB* p = n->child[i - 1];
            while (!p->leaf()) p = p->child[p->nr];
            n->key[i] = p->key[p->nr];
            k = n->key[i];
        }
        if (!removeFrom(n->child[i - 1], k)) return false;
        if (n->child[i - 1]->nr < MINK) fix(n, i - 1);
        return true;
    }

    // child c of n holds d - 1 keys: first try to borrow from the right sibling,
    // then from the left one; if both hold only d keys, the child is merged
    static void fix(NodeB* n, int c) {
        if (c < n->nr && n->child[c + 1]->nr > MINK) borrowRight(n, c);
        else if (c > 0 && n->child[c - 1]->nr > MINK) borrowLeft(n, c);
        else if (c < n->nr) merge(n, c);       // with the right sibling
        else merge(n, c - 1);                  // with the left sibling
    }

    // key[c + 1] moves down to the end of child c; the right sibling's first key
    // moves up
    static void borrowRight(NodeB* n, int c) {
        NodeB* a = n->child[c];
        NodeB* b = n->child[c + 1];
        a->nr++;
        a->key[a->nr] = n->key[c + 1];
        a->child[a->nr] = b->child[0];
        n->key[c + 1] = b->key[1];
        for (int j = 1; j < b->nr; j++) {
            b->key[j] = b->key[j + 1];
        }
        for (int j = 0; j < b->nr; j++) b->child[j] = b->child[j + 1];
        b->child[b->nr] = nullptr;
        b->nr--;
    }

    // key[c] moves down to the front of child c; the left sibling's last key moves up
    static void borrowLeft(NodeB* n, int c) {
        NodeB* a = n->child[c - 1];
        NodeB* b = n->child[c];
        for (int j = b->nr; j >= 1; j--) {
            b->key[j + 1] = b->key[j];
            b->child[j + 1] = b->child[j];
        }
        b->child[0] = a->child[a->nr];
        b->key[1] = n->key[c];
        b->nr++;
        n->key[c] = a->key[a->nr];
        a->child[a->nr] = nullptr;
        a->nr--;
    }

    // children c and c + 1 and the key key[c + 1] between them become one node
    static void merge(NodeB* n, int c) {
        NodeB* a = n->child[c];
        NodeB* b = n->child[c + 1];
        a->key[a->nr + 1] = n->key[c + 1];
        a->child[a->nr + 1] = b->child[0];
        for (int j = 1; j <= b->nr; j++) {     // the keys and children of b
            a->key[a->nr + 1 + j] = b->key[j];
            a->child[a->nr + 1 + j] = b->child[j];
        }
        a->nr += 1 + b->nr;
        for (int j = c + 1; j < n->nr; j++) {  // key[c + 1] and child[c + 1] leave n
            n->key[j] = n->key[j + 1];
            n->child[j] = n->child[j + 1];
        }
        n->child[n->nr] = nullptr;
        n->nr--;
        delete b;
    }

    static void inorder(const NodeB* n, std::vector<int>& out) {
        if (n == nullptr) return;
        for (int i = 1; i <= n->nr; i++) {
            inorder(n->child[i - 1], out);
            out.push_back(n->key[i]);
        }
        inorder(n->child[n->nr], out);
    }

    // goes down only into children that may hold keys in [a, b]
    static void inRange(const NodeB* n, int a, int b, std::vector<int>& out) {
        if (n == nullptr) return;
        int i = position(n, a);                // key[1..i-1] are smaller than a
        for (; i <= n->nr && n->key[i] <= b; i++) {
            inRange(n->child[i - 1], a, b, out);
            out.push_back(n->key[i]);
        }
        inRange(n->child[i - 1], a, b, out);   // the child after the last key <= b
    }

    static int countNodes(const NodeB* n) {
        if (n == nullptr) return 0;
        int c = 1;
        if (!n->leaf())
            for (int i = 0; i <= n->nr; i++) c += countNodes(n->child[i]);
        return c;
    }

    // the keys of n must lie in the interval (lo, hi); hasLo/hasHi say whether
    // the interval has a left and a right end
    static bool checkNode(const NodeB* n, int depth, bool isRoot, bool hasLo, int lo,
                          bool hasHi, int hi, int& leafLevel, int& total,
                          std::string& why) {
        if (n->nr < (isRoot ? 1 : MINK) || n->nr > MAXK) {
            why = "node " + str(n) + " holds " + std::to_string(n->nr) + " keys";
            return false;
        }
        for (int i = 1; i <= n->nr; i++) {
            if ((i > 1 && n->key[i] <= n->key[i - 1]) || (hasLo && n->key[i] <= lo) ||
                (hasHi && n->key[i] >= hi)) {
                why = "the keys of node " + str(n) + " are out of order";
                return false;
            }
        }
        total += n->nr;
        if (n->leaf()) {
            for (int i = 1; i <= n->nr; i++)
                if (n->child[i] != nullptr) {
                    why = "leaf " + str(n) + " has a child";
                    return false;
                }
            if (leafLevel < 0) leafLevel = depth;
            if (depth != leafLevel) {
                why = "the leaves are not all on the same level";
                return false;
            }
            return true;
        }
        for (int i = 0; i <= n->nr; i++) {
            if (n->child[i] == nullptr) {
                why = "internal node " + str(n) + " has a null child";
                return false;
            }
            bool l = i > 0 || hasLo, h = i < n->nr || hasHi;
            int a = i > 0 ? n->key[i] : lo, b = i < n->nr ? n->key[i + 1] : hi;
            if (!checkNode(n->child[i], depth + 1, false, l, a, h, b, leafLevel, total,
                           why))
                return false;
        }
        return true;
    }

    static void destroy(NodeB* n) {
        if (n == nullptr) return;
        for (int i = 0; i <= n->nr; i++) destroy(n->child[i]);
        delete n;
    }
};

// ---------------------------------------------------------------- menu

// reads a line; returns false at the end of the input
static bool readLine(const std::string& prompt, std::string& line) {
    std::cout << prompt;
    return static_cast<bool>(std::getline(std::cin, line));
}

static bool isSep(char c) {
    return c == ' ' || c == '\t' || c == ',' || c == ';' || c == '\r';
}

// reads the integers on a line; separators: spaces, commas, semicolons
static bool parseInts(const std::string& line, std::vector<int>& out) {
    out.clear();
    size_t i = 0;
    while (i < line.size()) {
        if (isSep(line[i])) { i++; continue; }
        size_t j = i;
        while (j < line.size() && !isSep(line[j])) j++;
        std::string tok = line.substr(i, j - i);
        errno = 0;
        char* end = nullptr;
        long v = std::strtol(tok.c_str(), &end, 10);
        if (*end != '\0' || errno == ERANGE || v < INT_MIN || v > INT_MAX) return false;
        out.push_back(static_cast<int>(v));
        i = j;
    }
    return true;
}

// reads exactly one integer
static bool readInt(const std::string& prompt, int& k, bool& eof) {
    std::string line;
    eof = !readLine(prompt, line);
    if (eof) return false;
    std::vector<int> v;
    if (!parseInts(line, v) || v.size() != 1) {
        std::cout << "Invalid value: enter a single integer.\n";
        return false;
    }
    k = v[0];
    return true;
}

// "1 key", "5 keys"
static std::string plural(int n, const char* one, const char* many) {
    return std::to_string(n) + " " + (n == 1 ? one : many);
}

static void printList(const std::vector<int>& v) {
    if (v.empty()) { std::cout << "  (no keys)\n"; return; }
    std::cout << " ";
    for (int k : v) std::cout << ' ' << k;
    std::cout << "\n  (" << plural(static_cast<int>(v.size()), "key", "keys")
              << ")\n";
}

int main() {
    BTree t;
    std::string line;
    std::cout << "B-trees - DSA lab\n"
              << "Order d = " << D << ": every node holds at most " << MAXK
              << " keys, and every node except the root holds at least " << MINK
              << ".\n";
    for (;;) {
        std::cout << "\n=== B-tree | " << plural(t.size(), "key", "keys") << " | "
                  << plural(t.nodes(), "node", "nodes") << " | height " << t.height()
                  << " ===\n"
                  << " 1. Insert a key\n"
                  << " 2. Insert several keys (on one line)\n"
                  << " 3. Search for a key\n"
                  << " 4. Delete a key\n"
                  << " 5. Print the tree level by level\n"
                  << " 6. Print the keys in increasing order\n"
                  << " 7. Print the keys in an interval [a, b]\n"
                  << " 8. Check the B-tree properties\n"
                  << " 9. Empty the tree\n"
                  << " 0. Exit\n";
        if (!readLine("Option: ", line)) break;
        std::vector<int> opt;
        if (!parseInts(line, opt) || opt.size() != 1 || opt[0] < 0 || opt[0] > 9) {
            std::cout << "Invalid option.\n";
            continue;
        }
        int o = opt[0];
        bool eof = false;
        int k = 0;
        if (o == 0) break;
        if (o == 1 || o == 3 || o == 4) {
            if (!readInt("Key: ", k, eof)) { if (eof) break; continue; }
            if (o == 1) {
                if (t.insert(k)) {
                    std::cout << "Inserted key " << k << ".\n";
                    t.print();
                } else {
                    std::cout << "Key " << k << " is already in the tree.\n";
                }
            } else if (o == 3) {
                std::vector<const NodeB*> path;
                bool found = t.search(k, &path);
                std::cout << "  path:";
                for (size_t i = 0; i < path.size(); i++)
                    std::cout << (i ? " -> " : " ") << BTree::str(path[i]);
                if (path.empty()) std::cout << " (the tree is empty)";
                std::cout << "\n  key " << k << (found ? " is" : " is not")
                          << " in the tree ("
                          << plural(static_cast<int>(path.size()), "node visited",
                                    "nodes visited")
                          << ")\n";
            } else {
                if (t.remove(k)) {
                    std::cout << "Deleted key " << k << ".\n";
                    t.print();
                } else {
                    std::cout << "Key " << k << " is not in the tree.\n";
                }
            }
        } else if (o == 2) {
            if (!readLine("Keys (separated by spaces or commas): ", line)) break;
            std::vector<int> v;
            if (!parseInts(line, v) || v.empty()) {
                std::cout << "Invalid list: enter integers.\n";
                continue;
            }
            int added = 0;
            for (int x : v) {
                if (t.insert(x)) added++;
                else std::cout << "  key " << x << " is already there\n";
            }
            std::cout << "Inserted " << plural(added, "key", "keys") << ".\n";
            t.print();
        } else if (o == 5) {
            t.print();
        } else if (o == 6) {
            printList(t.keys());
        } else if (o == 7) {
            int a = 0, b = 0;
            if (!readInt("Left end a: ", a, eof)) { if (eof) break; continue; }
            if (!readInt("Right end b: ", b, eof)) { if (eof) break; continue; }
            if (a > b) { int x = a; a = b; b = x; }
            printList(t.range(a, b));
        } else if (o == 8) {
            std::string why;
            if (t.check(why))
                std::cout << "The tree has all the properties of a B-tree.\n";
            else
                std::cout << "Property broken: " << why << ".\n";
        } else if (o == 9) {
            t.clear();
            std::cout << "The tree is now empty.\n";
        }
    }
    std::cout << "Goodbye!\n";
    return 0;
}
