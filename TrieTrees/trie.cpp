// trie_menu_en.cpp - Tries, Data Structures and Algorithms lab
//
// The variant with branch nodes and leaves:
//  - a branch node has 27 pointers: one for each of 'a'..'z' plus the blank pointer
//    (index 26), used by the word that ends exactly at that node;
//  - every word is kept whole in a leaf, placed as high in the tree as possible:
//    right below the first prefix that no other word shares.
//
// Build: g++ -std=c++11 -Wall -Wextra -o trie_menu trie_menu_en.cpp

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cctype>

const int ALPHABET = 26; // the letters 'a'..'z'
const int BLANK = ALPHABET; // index of the blank pointer
const int NR = ALPHABET + 1; // 27 pointers in a branch node

struct TrieNode {
    bool notLeaf; // true: branch node; false: leaf (holds a word)
    TrieNode *child[NR]; // children of a branch node; stay nullptr in leaves
    std::string word; // the word held by a leaf; empty in branch nodes

    explicit TrieNode(bool intern, const std::string &w = "")
        : notLeaf(intern), word(w) {
        for (int i = 0; i < NR; i++) child[i] = nullptr;
    }

    ~TrieNode() {
        for (int i = 0; i < NR; i++) delete child[i]; // deletes the whole subtree
    }
};

class Trie {
public:
    Trie() : root(new TrieNode(true)), nWords(0), nBranches(1) {
    }

    ~Trie() { delete root; }

    Trie(const Trie &) = delete; // the tree owns its nodes: no copies
    Trie &operator=(const Trie &) = delete;

    // true if s has only letters 'a'..'z' (an empty string is fine, as a prefix)
    static bool onlyLetters(const std::string &s) {
        for (char c: s)
            if (c < 'a' || c > 'z') return false;
        return true;
    }

    // a valid word is non-empty and holds only lowercase letters
    static bool isWord(const std::string &w) { return !w.empty() && onlyLetters(w); }

    // inserts w; returns false if w is invalid or already present
    bool insert(const std::string &w) {
        if (!isWord(w)) return false;
        TrieNode *node = root; // the current branch node
        size_t d = 0; // depth = letters consumed
        for (;;) {
            int i = slot(w, d);
            TrieNode *next = node->child[i];
            if (next == nullptr) {
                // free: the leaf goes here
                node->child[i] = new TrieNode(false, w);
                nWords++;
                return true;
            }
            if (next->notLeaf) {
                // branch node: go down
                node = next;
                d++;
                continue;
            }
            if (next->word == w) return false; // the word is already there
            // a leaf with another word: separate it from w with branch nodes
            TrieNode *old = next;
            do {
                TrieNode *branch = new TrieNode(true);
                node->child[i] = branch;
                nBranches++;
                node = branch;
                d++;
                i = slot(w, d);
            } while (i == slot(old->word, d)); // same letter: one more level
            node->child[slot(old->word, d)] = old;
            node->child[i] = new TrieNode(false, w);
            nWords++;
            return true;
        }
    }

    // searches for w; a leaf may come before the end of w: compare the whole word
    bool find(const std::string &w) const {
        if (!isWord(w)) return false;
        const TrieNode *node = root;
        for (size_t d = 0; ; d++) {
            const TrieNode *next = node->child[slot(w, d)];
            if (next == nullptr) return false; // no such branch
            if (!next->notLeaf) return next->word == w; // leaf: compare the word
            node = next; // branch node: go down
        }
    }

    // deletes w; returns false if w is not present
    bool remove(const std::string &w) {
        if (!isWord(w)) return false;
        std::vector<TrieNode *> path; // path[k]: node at depth k
        TrieNode *node = root;
        for (size_t d = 0; ; d++) {
            path.push_back(node);
            TrieNode *next = node->child[slot(w, d)];
            if (next == nullptr) return false;
            if (!next->notLeaf) {
                if (next->word != w) return false;
                node->child[slot(w, d)] = nullptr; // take the leaf out
                delete next;
                nWords--;
                break;
            }
            node = next;
        }
        // go up: a branch node left with only a leaf is replaced by that leaf
        for (size_t k = path.size() - 1; k > 0; k--) {
            TrieNode *branch = path[k];
            int only = onlyChild(branch);
            if (only < 0 || branch->child[only]->notLeaf) break;
            TrieNode *leaf = branch->child[only];
            branch->child[only] = nullptr; // detach the leaf
            path[k - 1]->child[slot(w, k - 1)] = leaf; // the leaf moves up
            delete branch;
            nBranches--;
        }
        return true;
    }

    // the words that start with p, in alphabetical order (empty p: all words)
    std::vector<std::string> wordsWithPrefix(const std::string &p) const {
        std::vector<std::string> out;
        if (!onlyLetters(p)) return out;
        const TrieNode *node = root;
        for (size_t d = 0; d < p.size(); d++) {
            const TrieNode *next = node->child[p[d] - 'a'];
            if (next == nullptr) return out; // no word has this prefix
            if (!next->notLeaf) {
                // a leaf on the path of p
                if (next->word.compare(0, p.size(), p) == 0) out.push_back(next->word);
                return out;
            }
            node = next;
        }
        collect(node, out); // the subtree has prefix p
        return out;
    }

    // draws the tree; '_' marks the blank pointer
    void print() const {
        std::cout << "(root)\n";
        if (nWords == 0) std::cout << "    the tree is empty\n";
        else draw(root, "");
    }

    int words() const { return nWords; }
    int branches() const { return nBranches; }

    void clear() {
        delete root;
        root = new TrieNode(true);
        nWords = 0;
        nBranches = 1;
    }

private:
    TrieNode *root; // the root is always a branch node
    int nWords; // number of words = number of leaves
    int nBranches; // number of branch nodes, the root included

    // the child followed at depth d: letter w[d], or blank once w has ended
    static int slot(const std::string &w, size_t d) {
        return d < w.size() ? w[d] - 'a' : BLANK;
    }

    // children are visited blank first, then 'a'..'z': alphabetical order
    static int nth(int k) { return k == 0 ? BLANK : k - 1; }

    // index of the only child of n, or -1 if n has several children
    static int onlyChild(const TrieNode *n) {
        int found = -1;
        for (int i = 0; i < NR; i++)
            if (n->child[i] != nullptr) {
                if (found != -1) return -1;
                found = i;
            }
        return found;
    }

    // preorder traversal: appends every word in the subtree of n to out
    void collect(const TrieNode *n, std::vector<std::string> &out) const {
        if (!n->notLeaf) {
            out.push_back(n->word);
            return;
        }
        for (int k = 0; k < NR; k++)
            if (n->child[nth(k)] != nullptr) collect(n->child[nth(k)], out);
    }

    void draw(const TrieNode *n, const std::string &indent) const {
        std::vector<int> kids;
        for (int k = 0; k < NR; k++)
            if (n->child[nth(k)] != nullptr) kids.push_back(nth(k));
        for (size_t j = 0; j < kids.size(); j++) {
            bool last = (j + 1 == kids.size());
            const TrieNode *c = n->child[kids[j]];
            char label = (kids[j] == BLANK) ? '_' : char('a' + kids[j]);
            std::cout << indent << (last ? "`-- " : "+-- ") << label;
            if (!c->notLeaf) std::cout << " : " << c->word << '\n';
            else {
                std::cout << '\n';
                draw(c, indent + (last ? "    " : "|   "));
            }
        }
    }
};

// ---------------------------------------------------------------- menu

// reads a line; returns false at the end of the input
static bool readLine(const std::string &prompt, std::string &line) {
    std::cout << prompt;
    return static_cast<bool>(std::getline(std::cin, line));
}

// "1 word", "5 words"
static std::string plural(int n, const std::string &one, const std::string &many) {
    return std::to_string(n) + ' ' + (n == 1 ? one : many);
}

// lowercase letters, no surrounding spaces; capital letters are accepted
static std::string normalize(const std::string &s) {
    std::string r;
    for (char c: s)
        if (!std::isspace(static_cast<unsigned char>(c)))
            r += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return r;
}

static void printList(const std::vector<std::string> &v) {
    if (v.empty()) {
        std::cout << "  (no words)\n";
        return;
    }
    std::cout << " ";
    for (const std::string &w: v) std::cout << ' ' << w;
    int n = static_cast<int>(v.size());
    std::cout << "\n  (" << plural(n, "word", "words") << ")\n";
}

int main() {
    Trie t;
    std::string line;
    std::cout << "Tries - DSA lab\n"
            "Words use only the letters a-z "
            "(capital letters are turned into lowercase).\n";
    for (;;) {
        std::cout << "\n=== trie | " << plural(t.words(), "word", "words") << " | "
                << plural(t.branches(), "branch node", "branch nodes")
                << " ===\n"
                << " 1. Insert a word\n"
                << " 2. Insert several words (on one line)\n"
                << " 3. Search for a word\n"
                << " 4. Delete a word\n"
                << " 5. List the words in alphabetical order\n"
                << " 6. List the words that start with a prefix\n"
                << " 7. Draw the tree\n"
                << " 8. Clear the tree\n"
                << " 0. Exit\n";
        if (!readLine("Option: ", line)) break;
        std::string opt = normalize(line);
        if (opt == "0") break;

        if (opt == "1" || opt == "3" || opt == "4") {
            if (!readLine("Word: ", line)) break;
            std::string w = normalize(line);
            if (!Trie::isWord(w)) {
                std::cout << "Invalid word: use only the letters a-z.\n";
                continue;
            }
            if (opt == "1") {
                if (t.insert(w)) {
                    std::cout << "Inserted \"" << w << "\".\n";
                    t.print();
                } else {
                    std::cout << "\"" << w << "\" is already in the tree.\n";
                }
            } else if (opt == "3") {
                std::cout << (t.find(w) ? "Found: \"" : "Not found: \"") << w << "\"\n";
            } else {
                if (t.remove(w)) {
                    std::cout << "Deleted \"" << w << "\".\n";
                    t.print();
                } else {
                    std::cout << "\"" << w << "\" is not in the tree.\n";
                }
            }
        } else if (opt == "2") {
            const char *ask = "Words separated by spaces (e.g. car cart care cat): ";
            if (!readLine(ask, line)) break;
            for (char &c: line) if (c == ',' || c == ';') c = ' ';
            std::istringstream in(line);
            std::string token;
            int added = 0;
            while (in >> token) {
                std::string w = normalize(token);
                if (!Trie::isWord(w))
                    std::cout << "  skipped (invalid): " << token << '\n';
                else if (t.insert(w))
                    added++;
                else
                    std::cout << "  already there: " << w << '\n';
            }
            std::cout << "Inserted " << plural(added, "word", "words") << ".\n";
            t.print();
        } else if (opt == "5") {
            printList(t.wordsWithPrefix(""));
        } else if (opt == "6") {
            if (!readLine("Prefix: ", line)) break;
            std::string p = normalize(line);
            if (!Trie::onlyLetters(p)) {
                std::cout << "Invalid prefix: use only the letters a-z.\n";
                continue;
            }
            printList(t.wordsWithPrefix(p));
        } else if (opt == "7") {
            t.print();
        } else if (opt == "8") {
            t.clear();
            std::cout << "The tree is now empty.\n";
        } else {
            std::cout << "Unknown option. Choose a number from 0 to 8.\n";
        }
    }
    return 0;
}
