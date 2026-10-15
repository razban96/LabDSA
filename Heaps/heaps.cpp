// heap_menu_en.cpp — menu-driven binary heap (min-heap or max-heap)
// Data Structures and Algorithms lab, University of Craiova
//
// Starts from the MinHeap class on GeeksforGeeks (harr, capacity, heap_size,
// parent/left/right, insertKey, decreaseKey, extractMin, deleteKey, MinHeapify),
// https://www.geeksforgeeks.org/dsa/binary-heap/, and adds:
//  - a mirrored max-heap version: every min... method has a max... twin
//  - an interactive menu, with the heap type chosen at start-up and switchable later
//  - a console drawing of the tree and a "step by step" mode that shows
//    every swap made by sift-up and heapify
// The elements live in a std::vector, so the capacity is no longer fixed.
//
// Build: g++ -std=c++11 -Wall -Wextra -o heap_menu_en heap_menu_en.cpp

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

const int DRAW_LEVELS = 6; // taller trees are printed level by level

enum class HeapType { Min, Max };

class Heap {
public:
    explicit Heap(HeapType t) : type(t) {}
    HeapType getType() const { return type; }
    int size() const { return static_cast<int>(harr.size()); }
    int at(int i) const { return harr[i]; }
    bool getTrace() const { return trace; }
    void setTrace(bool on) { trace = on; }

    void insertKey(int k);
    bool getRoot(int& out) const;
    bool extractMin(int& out);
    bool extractMax(int& out);
    bool extractRoot(int& out);
    bool decreaseKey(int i, int newVal);
    bool increaseKey(int i, int newVal);
    bool changeKey(int i, int newVal);
    bool deleteKey(int i, int& out);
    void buildHeap();
    void load(const std::vector<int>& v);
    void switchType();
    bool isHeap() const;
    std::vector<int> sorted() const;

    void printArray() const;
    void printTree() const;

private:
    std::vector<int> harr; // the heap elements
    HeapType type;        // HeapType::Min or HeapType::Max
    bool trace = false;   // true: print every swap

    static int parent(int i) { return (i - 1) / 2; }
    static int left(int i)   { return 2 * i + 1; }
    static int right(int i)  { return 2 * i + 2; }

    void minSiftUp(int i);
    void maxSiftUp(int i);
    void minHeapify(int i, int n);
    void maxHeapify(int i, int n);
    void siftUp(int i);
    void heapify(int i, int n);
    void traceSwap(int i, int j) const;
};

//
// Helpers
//
static const char* typeName(HeapType t) {
    return t == HeapType::Min ? "min-heap" : "max-heap";
}

static const char* rootName(HeapType t) {
    return t == HeapType::Min ? "minimum" : "maximum";
}

// Number of characters x takes when printed (including the minus sign).
static int numWidth(int x) {
    return static_cast<int>(std::to_string(x).size());
}

// Called right after a swap between positions i and j.
void Heap::traceSwap(int i, int j) const {
    if (!trace) return;
    std::cout << "  swap h[" << i << "] = " << harr[j]
              << " with h[" << j << "] = " << harr[i] << '\n';
    printTree();
}

//
// Moving up (sift-up): the min / max pair
//
void Heap::minSiftUp(int i) {
    while (i > 0 && harr[i] < harr[parent(i)]) {
        std::swap(harr[i], harr[parent(i)]);
        traceSwap(i, parent(i));
        i = parent(i);
    }
}

void Heap::maxSiftUp(int i) {
    while (i > 0 && harr[i] > harr[parent(i)]) {
        std::swap(harr[i], harr[parent(i)]);
        traceSwap(i, parent(i));
        i = parent(i);
    }
}

//
// Moving down (heapify / sift-down): the min / max pair
// Works on the first n elements and assumes the subtrees of i
// are already heaps.
//
void Heap::minHeapify(int i, int n) {
    int l = left(i), r = right(i), smallest = i;
    if (l < n && harr[l] < harr[smallest]) smallest = l;
    if (r < n && harr[r] < harr[smallest]) smallest = r;
    if (smallest != i) {
        std::swap(harr[i], harr[smallest]);
        traceSwap(i, smallest);
        minHeapify(smallest, n);
    }
}

void Heap::maxHeapify(int i, int n) {
    int l = left(i), r = right(i), largest = i;
    if (l < n && harr[l] > harr[largest]) largest = l;
    if (r < n && harr[r] > harr[largest]) largest = r;
    if (largest != i) {
        std::swap(harr[i], harr[largest]);
        traceSwap(i, largest);
        maxHeapify(largest, n);
    }
}

// Pick the version that matches the heap type.
void Heap::siftUp(int i) {
    if (type == HeapType::Min) minSiftUp(i); else maxSiftUp(i);
}

void Heap::heapify(int i, int n) {
    if (type == HeapType::Min) minHeapify(i, n); else maxHeapify(i, n);
}

//
// Operations
//
void Heap::insertKey(int k) {
    harr.push_back(k); // first free position; the vector grows by itself
    siftUp(size() - 1);
}

// The root is the minimum (min-heap) or the maximum (max-heap).
bool Heap::getRoot(int& out) const {
    if (harr.empty()) return false;
    out = harr[0];
    return true;
}

bool Heap::extractMin(int& out) {
    if (harr.empty()) return false;
    out = harr[0];
    harr[0] = harr.back(); // the last element moves to the root
    harr.pop_back();
    minHeapify(0, size());
    return true;
}

bool Heap::extractMax(int& out) {
    if (harr.empty()) return false;
    out = harr[0];
    harr[0] = harr.back();
    harr.pop_back();
    maxHeapify(0, size());
    return true;
}

bool Heap::extractRoot(int& out) {
    return type == HeapType::Min ? extractMin(out) : extractMax(out);
}

// Min-heap: the key may only decrease, so the element may only move up.
bool Heap::decreaseKey(int i, int newVal) {
    if (i < 0 || i >= size() || newVal > harr[i]) return false;
    harr[i] = newVal;
    minSiftUp(i);
    return true;
}

// Max-heap: the key may only increase, so the element may only move up.
bool Heap::increaseKey(int i, int newVal) {
    if (i < 0 || i >= size() || newVal < harr[i]) return false;
    harr[i] = newVal;
    maxSiftUp(i);
    return true;
}

// Change in either direction: towards the root with decreaseKey / increaseKey,
// otherwise the element moves down with heapify.
bool Heap::changeKey(int i, int newVal) {
    if (i < 0 || i >= size()) return false;
    if (type == HeapType::Min) {
        if (newVal <= harr[i]) return decreaseKey(i, newVal);
        harr[i] = newVal;
        minHeapify(i, size());
    } else {
        if (newVal >= harr[i]) return increaseKey(i, newVal);
        harr[i] = newVal;
        maxHeapify(i, size());
    }
    return true;
}

// Delete the element at position i: replace it with the last element
// and restore the order in whichever direction is needed.
// Classic alternative: decreaseKey(i, INT_MIN) + extractMin for a min-heap,
// or increaseKey(i, INT_MAX) + extractMax for a max-heap.
bool Heap::deleteKey(int i, int& out) {
    if (i < 0 || i >= size()) return false;
    out = harr[i];
    int last = harr.back();
    harr.pop_back();
    if (i < size()) changeKey(i, last);
    return true;
}

// Build-heap in O(n): heapify the internal nodes, from the last one up to the root.
void Heap::buildHeap() {
    for (int i = size() / 2 - 1; i >= 0; i--) {
        heapify(i, size());
    }
}

void Heap::load(const std::vector<int>& v) {
    harr = v;
    buildHeap();
}

void Heap::switchType() {
    type = (type == HeapType::Min) ? HeapType::Max : HeapType::Min;
    buildHeap();
}

bool Heap::isHeap() const {
    for (int i = 1; i < size(); i++) {
        int p = harr[parent(i)], c = harr[i];
        if (type == HeapType::Min ? p > c : p < c) return false;
    }
    return true;
}

// In-place heapsort on a copy, so the heap in the menu stays unchanged.
// min-heap -> descending order, max-heap -> ascending order.
std::vector<int> Heap::sorted() const {
    Heap copy = *this;
    copy.trace = false;
    for (int end = copy.size() - 1; end > 0; end--) {
        std::swap(copy.harr[0], copy.harr[end]);
        copy.heapify(0, end); // the root goes to its final place
    }
    return copy.harr;
}

//
// Output
//
void Heap::printArray() const {
    if (harr.empty()) { std::cout << "The heap is empty.\n"; return; }
    int w = numWidth(size() - 1);
    for (int x : harr) w = std::max(w, numWidth(x));
    w++;
    std::cout << "  index:";
    for (int i = 0; i < size(); i++) std::cout << std::setw(w) << i;
    std::cout << "\n  value:";
    for (int x : harr) std::cout << std::setw(w) << x;
    std::cout << '\n';
}

// Print the line without its trailing spaces.
static void printTrimmed(const std::string& s) {
    std::size_t end = s.find_last_not_of(' ');
    std::cout << "  " << (end == std::string::npos ? "" : s.substr(0, end + 1)) << '\n';
}

// Tree drawing: each node is centred in its own "slice" of the width.
// Level d has 2^d slices; a leaf on the last level gets cw characters.
void Heap::printTree() const {
    int n = size(), height = 0, cw = 4;
    if (n == 0) { std::cout << "  (empty heap)\n"; return; }
    while ((2 << height) - 1 < n) height++; // height = floor(log n)

    for (int x : harr) cw = std::max(cw, numWidth(x) + 2);

    if (height >= DRAW_LEVELS) { // too wide for a console
        for (int d = 0, first = 0; first < n; d++, first = 2 * first + 1) {
            std::cout << "  level " << d << ":";
            for (int i = first; i < n && i <= 2 * first; i++)
                std::cout << " " << harr[i];
            std::cout << '\n';
        }
        return;
    }

    int total = (1 << height) * cw;
    for (int d = 0; d <= height; d++) {
        int slots = 1 << d, sw = total / slots, first = slots - 1;
        std::string line(total, ' ');
        for (int k = 0; k < slots && first + k < n; k++) {
            std::string s = std::to_string(harr[first + k]);
            int len = static_cast<int>(s.size());
            line.replace(k * sw + (sw - len) / 2, s.size(), s);
        }
        printTrimmed(line);

        if (d < height) {
            std::string edges(total, ' ');
            // edges towards the children
            for (int k = 0; k < slots && first + k < n; k++) {
                int i = first + k;
                if (left(i) < n)  edges[k * sw + 3 * sw / 8] = '/';
                if (right(i) < n) edges[k * sw + 5 * sw / 8] = '\\';
            }
            printTrimmed(edges);
        }
    }
}

//
// Keyboard input
//
static std::ifstream testInput;

static bool readLine(const std::string& prompt, std::string& line) {
    std::cout << prompt << std::flush;
    while (std::getline(testInput, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        return true;
    }
    return false; // EOF
}

// Repeats until it gets a valid integer. Returns false on EOF.
static bool readInt(const std::string& prompt, int& out) {
    std::string line;
    while (readLine(prompt, line)) {
        std::istringstream in(line);
        int v;
        char extra;
        if (in >> v && !(in >> extra)) { out = v; return true; }
        std::cout << "  Invalid value, please enter an integer.\n";
    }
    return false;
}

// Reads numbers separated by spaces or commas from one line.
// Returns how many it read, or -1 on EOF.
static int readList(const std::string& prompt, std::vector<int>& v) {
    std::string line;
    if (!readLine(prompt, line)) return -1;
    for (char& c : line)
        if (c == ',' || c == ';') c = ' ';
    std::istringstream in(line);
    v.clear();
    int x;
    for (;;) {
        if (in >> x) v.push_back(x);
        else if (in.eof()) break;
        else { in.clear(); in.ignore(1); } // skip one invalid character
    }
    return static_cast<int>(v.size());
}

//
// Menu
//
static void printMenu(const Heap& h) {
    bool mn = (h.getType() == HeapType::Min);
    std::cout << "\n=== " << typeName(h.getType()) << " (" << h.size()
              << " elements) | step-by-step: " << (h.getTrace() ? "ON" : "OFF")
              << " ===\n"
              << "  1. Insert a value\n"
              << "  2. Show the root (" << rootName(h.getType()) << ")\n"
              << "  3. Extract the root (" << (mn ? "extractMin" : "extractMax") << ")\n"
              << "  4. " << (mn ? "Decrease" : "Increase") << " the key at an index ("
              << (mn ? "decreaseKey" : "increaseKey") << ")\n"
              << "  5. Change the key at an index, in either direction\n"
              << "  6. Delete the element at an index\n"
              << "  7. Build a heap from a list of values (buildHeap)\n"
              << "  8. Heapsort (" << (mn ? "descending" : "ascending") << " order)\n"
              << "  9. Show the array\n"
              << " 10. Show the tree\n"
              << " 11. Check the heap property\n"
              << " 12. Switch to " << (mn ? "max-heap" : "min-heap")
              << " (rebuilds the heap)\n"
              << " 13. Turn step-by-step mode " << (h.getTrace() ? "off" : "on") << "\n"
              << "  0. Exit\n";
}

static bool validIndex(const Heap& h, int i) {
    if (i >= 0 && i < h.size()) return true;
    if (h.size() == 0)
        std::cout << "The heap is empty.\n";
    else
        std::cout << "Invalid index: it must be between 0 and "
                  << h.size() - 1 << ".\n";
    return false;
}

int main() {
    int choice, t, i, x;
    std::vector<int> v;
    bool running = true;

    testInput.open("input.txt");
    if (!testInput) {
        std::cerr << "Could not open input.txt. Run the program from the folder "
                     "that contains heaps.cpp and input.txt.\n";
        return 1;
    }

    std::cout << "Binary heap - DSA lab\n"
              << "  1. min-heap (parent <= children, the root is the minimum)\n"
              << "  2. max-heap (parent >= children, the root is the maximum)\n";

    do {
        if (!readInt("Choose the type (1 or 2): ", t)) return 0;
    } while (t != 1 && t != 2);

    Heap h(t == 1 ? HeapType::Min : HeapType::Max);

    while (running) {
        printMenu(h);
        if (!readInt("Option: ", choice)) break;

        switch (choice) {
            case 0:
                running = false;
                break;
            case 1:
                if (!readInt("Value: ", x)) { running = false; break; }
                h.insertKey(x);
                std::cout << "Inserted " << x << ". Result:\n";
                h.printTree();
                break;
            case 2:
                if (h.getRoot(x))
                    std::cout << "The root (" << rootName(h.getType()) << ") is "
                              << x << ".\n";
                else
                    std::cout << "The heap is empty.\n";
                break;
            case 3:
                if (h.extractRoot(x)) {
                    std::cout << "Extracted " << x << ". Result:\n";
                    h.printTree();
                } else {
                    std::cout << "The heap is empty.\n";
                }
                break;
            case 4: {
                if (!readInt("Index: ", i)) { running = false; break; }
                if (!validIndex(h, i)) break;
                if (!readInt("New value: ", x)) { running = false; break; }
                bool mn = (h.getType() == HeapType::Min);
                if (mn ? h.decreaseKey(i, x) : h.increaseKey(i, x)) {
                    std::cout << "Result:\n";
                    h.printTree();
                } else {
                    std::cout << "In a " << typeName(h.getType()) << " "
                              << (mn ? "decreaseKey" : "increaseKey") << " needs a value "
                              << (mn ? "smaller" : "greater") << " than or equal to "
                              << h.at(i) << "\nUse option 5 instead.\n";
                }
                break;
            }
            case 5:
                if (!readInt("Index: ", i)) { running = false; break; }
                if (!validIndex(h, i)) break;
                if (!readInt("New value: ", x)) { running = false; break; }
                h.changeKey(i, x);
                std::cout << "Result:\n";
                h.printTree();
                break;
            case 6:
                if (!readInt("Index: ", i)) { running = false; break; }
                if (!validIndex(h, i)) break;
                h.deleteKey(i, x);
                std::cout << "Deleted " << x << ". Result:\n";
                h.printTree();
                break;
            case 7: {
                int n = readList("Values (e.g. 41 7 23 15 3 36 9 28 12 50): ", v);
                if (n < 0) { running = false; break; }
                if (n == 0) { std::cout << "No values were read.\n"; break; }
                h.load(v);
                std::cout << "Built a " << typeName(h.getType()) << " with " << n
                          << " elements:\n";
                h.printArray();
                h.printTree();
                break;
            }
            case 8:
                if (h.size() == 0) { std::cout << "The heap is empty.\n"; break; }
                std::cout << "Sorted "
                          << (h.getType() == HeapType::Min ? "descending" : "ascending")
                          << ": ";
                for (int y : h.sorted()) std::cout << y << " ";
                std::cout << '\n';
                break;
            case 9:
                h.printArray();
                break;
            case 10:
                h.printTree();
                break;
            case 11:
                std::cout << "The array "
                          << (h.isHeap() ? "satisfies" : "does NOT satisfy")
                          << " the " << typeName(h.getType()) << " property.\n";
                break;
            case 12:
                h.switchType();
                std::cout << "Switched to " << typeName(h.getType()) << ". Result:\n";
                h.printTree();
                break;
            case 13:
                h.setTrace(!h.getTrace());
                std::cout << "Step-by-step mode is "
                          << (h.getTrace() ? "on" : "off") << ".\n";
                break;
            default:
                std::cout << "Unknown option.\n";
        }
    }

    std::cout << "Goodbye!\n";
    return 0;
}