// Single-source shortest paths (SSSP) - DSA lab
//
// For a directed graph with arc costs and a source vertex s, the program computes the
// shortest distance d[v] from s to every vertex v and the predecessor of v on a
// shortest path:
// - Dijkstra without a heap: the minimum is found by a linear scan of d, O(n^2 + m);
// - Dijkstra with a min-heap: the reached, unfinished vertices sit in a min-heap,
//   O((n + m) log n);
// - Bellman-Ford: at most n - 1 passes over all the arcs, O(n m); accepts negative
//   costs and detects the negative-cost cycles reachable from s.
// Dijkstra gives correct distances only if no arc has a negative cost.
//
// A graph file: n and m on the first line, then m lines "u v c", the arc u -> v of cost
// c. Vertices are numbered from 1 to n; the numbering 0..n-1 is accepted as well. Lines
// that start with # are comments.
//
// Build: g++ -std=c++11 -O2 -Wall -Wextra -o sssp sssp_en.cpp
// Run: ./sssp or ./sssp graf_mic.txt

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

typedef long long Dist;
const Dist INF = LLONG_MAX / 4;
const Dist NEG_INF = -INF;
const int MAX_COST = 1000000000;

struct Arc { int to, cost; };
struct Edge { int from, to, cost; };

struct Graph {
    int n = 0;
    int base = 1;
    std::string note;
    std::vector<Edge> edges;
    std::vector<std::vector<Arc>> adj;

    void reset(int vertices, int numbering) {
        n = vertices;
        base = numbering;
        note.clear();
        edges.clear();
        adj.assign(n, std::vector<Arc>());
    }

    void add(int u, int v, int c) {
        edges.push_back(Edge{u, v, c});
        adj[u].push_back(Arc{v, c});
    }

    int m() const { return static_cast<int>(edges.size()); }
    int label(int v) const { return v + base; }

    bool hasNegative() const {
        for (const Edge& e : edges)
            if (e.cost < 0) return true;
        return false;
    }
};

struct Result {
    std::vector<Dist> dist;
    std::vector<int> pred;
    std::vector<int> order;
    std::vector<int> cycle;
    long long relax = 0;
    long long updates = 0;
    int passes = 0;
    bool negCycle = false;

    explicit Result(int n) : dist(n, INF), pred(n, -1) {}
};

// ---------------------------------------------------------------- min-heap

// Min-heap of vertices ordered by (d[v], v): on equal distances the vertex with the
// smaller index comes out first, as in the version without a heap. pos[v] is the
// position of v in the heap or -1, so that decrease finds the vertex in O(1). The
// children of position i are 2i + 1 and 2i + 2.
class MinHeap {
public:
    MinHeap(int n, const std::vector<Dist>& d) : d_(d), pos_(n, -1) { h_.reserve(n); }
    bool empty() const { return h_.empty(); }
    bool contains(int v) const { return pos_[v] >= 0; }

    void push(int v) { // inserts v, with key d[v]
        h_.push_back(v);
        up(static_cast<int>(h_.size()) - 1);
    }

    void decrease(int v) { up(pos_[v]); } // d[v] has just dropped: v sifts up

    int pop() { // extracts the vertex with the smallest distance
        int top = h_, last = h_.back();
        h_.pop_back();
        pos_[top] = -1;
        if (!h_.empty()) {
            h_ = last;
            down(0);
        }
        return top;
    }

private:
    const std::vector<Dist>& d_;
    std::vector<int> h_, pos_;

    bool less(int a, int b) const { return d_[a] < d_[b] || (d_[a] == d_[b] && a < b); }
    void place(int v, int i) { h_[i] = v; pos_[v] = i; }

    void up(int i) {
        int v = h_[i];
        while (i > 0 && less(v, h_[(i - 1) / 2])) {
            place(h_[(i - 1) / 2], i);
            i = (i - 1) / 2;
        }
        place(v, i);
    }

    void down(int i) {
        int v = h_[i], size = static_cast<int>(h_.size());
        for (;;) {
            int c = 2 * i + 1;
            if (c >= size) break;
            if (c + 1 < size && less(h_[c + 1], h_[c])) c++;
            if (!less(h_[c], v)) break;
            place(h_[c], i);
            i = c;
        }
        place(v, i);
    }
};

// ---------------------------------------------------------------- algorithms

// Dijkstra without a heap: at every step a linear scan picks the unfinished vertex with
// the smallest distance; then the arcs that leave it are relaxed.
Result dijkstraArray(const Graph& g, int s) {
    Result r(g.n);
    std::vector<char> done(g.n, 0);
    r.dist[s] = 0;
    for (int step = 0; step < g.n; step++) {
        int u = -1;
        for (int v = 0; v < g.n; v++)
            if (!done[v] && r.dist[v] < INF && (u < 0 || r.dist[v] < r.dist[u])) u = v;
        if (u < 0) break; // the other vertices are unreachable
        done[u] = 1;
        r.order.push_back(u);
        for (const Arc& a : g.adj[u]) {
            r.relax++;
            if (!done[a.to] && r.dist[u] + a.cost < r.dist[a.to]) {
                r.dist[a.to] = r.dist[u] + a.cost;
                r.pred[a.to] = u;
                r.updates++;
            }
        }
    }
    return r;
}

// Dijkstra with a min-heap: the minimum is extracted from the heap in O(log n); a
// vertex enters the heap when it is first reached and sifts up when its distance
// drops.
Result dijkstraHeap(const Graph& g, int s) {
    Result r(g.n);
    std::vector<char> done(g.n, 0);
    MinHeap h(g.n, r.dist);
    r.dist[s] = 0;
    h.push(s);
    while (!h.empty()) {
        int u = h.pop();
        done[u] = 1;
        r.order.push_back(u);
        for (const Arc& a : g.adj[u]) {
            r.relax++;
            if (!done[a.to] && r.dist[u] + a.cost < r.dist[a.to]) {
                r.dist[a.to] = r.dist[u] + a.cost;
                r.pred[a.to] = u;
                r.updates++;
                if (h.contains(a.to)) h.decrease(a.to);
                else h.push(a.to);
            }
        }
    }
    return r;
}

// the vertices reachable from a vertex affected by a negative cycle have distance -inf
static void markNegInf(const Graph& g, std::vector<int> stack, std::vector<Dist>& dist) {
    for (int v : stack) dist[v] = NEG_INF;
    while (!stack.empty()) {
        int u = stack.back();
        stack.pop_back();
        for (const Arc& a : g.adj[u])
            if (dist[a.to] != NEG_INF) {
                dist[a.to] = NEG_INF;
                stack.push_back(a.to);
            }
    }
}

// Bellman-Ford: every pass relaxes all the arcs, in the order of the list. After pass
// k, d[v] is at most the cost of any path from s to v with at most k arcs, so n - 1
// passes are enough. If a pass changes nothing, the distances are final.
Result bellmanFord(const Graph& g, int s) {
    Result r(g.n);
    r.dist[s] = 0;
    bool stable = false;
    for (int pass = 1; pass < g.n && !stable; pass++) {
        stable = true;
        for (const Edge& e : g.edges) {
            r.relax++;
            if (r.dist[e.from] != INF && r.dist[e.from] + e.cost < r.dist[e.to]) {
                r.dist[e.to] = r.dist[e.from] + e.cost;
                r.pred[e.to] = e.from;
                r.updates++;
                stable = false;
            }
        }
        r.passes = pass;
    }
    if (stable) return r;
    // pass n: an arc that can still be relaxed shows a negative cycle
    std::vector<int> bad;
    for (const Edge& e : g.edges) {
        r.relax++;
        if (r.dist[e.from] != INF && r.dist[e.from] + e.cost < r.dist[e.to]) {
            r.dist[e.to] = r.dist[e.from] + e.cost;
            r.pred[e.to] = e.from;
            bad.push_back(e.to);
        }
    }
    if (bad.empty()) return r;
    r.negCycle = true;
    int x = bad.back(); // after n steps back, x lies on the cycle
    for (int i = 0; i < g.n && x >= 0; i++) x = r.pred[x];
    if (x >= 0) {
        int v = x;
        do {
            r.cycle.push_back(v);
            v = r.pred[v];
        } while (v != x && v >= 0 && r.cycle.size() <= static_cast<size_t>(g.n));
        std::reverse(r.cycle.begin(), r.cycle.end()); // in the direction of the arcs
        std::rotate(r.cycle.begin(),                  // start at the smallest vertex
                    std::min_element(r.cycle.begin(), r.cycle.end()), r.cycle.end());
        r.cycle.push_back(r.cycle);
    }
    markNegInf(g, bad, r.dist);
    return r;
}

// the shortest path from s to t, rebuilt from the predecessors; empty if none is finite
static std::vector<int> pathTo(const Result& r, int t) {
    std::vector<int> p;
    if (r.dist[t] == INF || r.dist[t] == NEG_INF) return p;
    for (int v = t; v != -1; v = r.pred[v]) p.push_back(v);
    std::reverse(p.begin(), p.end());
    return p;
}

// ---------------------------------------------------------------- files

// "1 arc", "5 arcs"
static std::string plural(long long n, const char* one, const char* many) {
    return std::to_string(n) + " " + (n == 1 ? one : many);
}

static bool isSep(char c) {
    return c == ' ' || c == '\t' || c == ',' || c == ';' || c == '\r';
}

// reads the integers of a line; separators: spaces, commas, semicolons
static bool parseNums(const std::string& line, std::vector<long long>& out) {
    out.clear();
    size_t i = 0;
    while (i < line.size()) {
        if (isSep(line[i])) { i++; continue; }
        size_t j = i;
        while (j < line.size() && !isSep(line[j])) j++;
        std::string tok = line.substr(i, j - i);
        errno = 0;
        char* end = nullptr;
        long long v = std::strtoll(tok.c_str(), &end, 10);
        if (*end != '\0' || errno == ERANGE) return false;
        out.push_back(v);
        i = j;
    }
    return true;
}

static bool blankOrComment(const std::string& line) {
    for (char c : line) {
        if (c == '#') return true;
        if (!isSep(c)) return false;
    }
    return true;
}

// reads the graph from a file; on error returns false and explains in err
static bool loadGraph(const std::string& file, Graph& g, std::string& err) {
    std::ifstream in(file.c_str());
    if (!in) { err = "the file cannot be opened"; return false; }
    std::string line;
    std::vector<long long> v;
    long long n = -1, m = -1;
    std::vector<Edge> list;
    std::vector<int> lineOf;
    int no = 0, lo = INT_MAX;
    while (std::getline(in, line)) {
        no++;
        if (blankOrComment(line)) continue;
        std::string where = "line " + std::to_string(no) + ": ";
        if (!parseNums(line, v)) {
            err = where + "integers were expected";
            return false;
        }
        if (n < 0) {
            if (v.size() != 2) { err = where + "n and m were expected"; return false; }
            n = v;
            m = v;
            if (n < 1 || n > 1000000) {
                err = where + "n must be between 1 and 10^6";
                return false;
            }
            if (m < 0 || m > 50000000) {
                err = where + "m must be between 0 and 5*10^7";
                return false;
            }
            list.reserve(static_cast<size_t>(m));
            continue;
        }
        if (v.size() != 3) { err = where + "an arc was expected: u v c"; return false; }
        if (static_cast<long long>(list.size()) == m) {
            err = where + "the file has more than m = " + std::to_string(m) + " arcs";
            return false;
        }
        if (v < INT_MIN || v > INT_MAX || v < INT_MIN || v > INT_MAX) {
            err = where + "vertex out of range";
            return false;
        }
        if (v < -MAX_COST || v > MAX_COST) {
            err = where + "the cost must be between -10^9 and 10^9";
            return false;
        }
        Edge e = {static_cast<int>(v), static_cast<int>(v), static_cast<int>(v)};
        lo = std::min(lo, std::min(e.from, e.to));
        list.push_back(e);
        lineOf.push_back(no);
    }
    if (n < 0) { err = "the file does not hold n and m"; return false; }
    if (static_cast<long long>(list.size()) != m) {
        long long k = static_cast<long long>(list.size());
        err = "the file has " + plural(k, "arc", "arcs") + ", not m = " + std::to_string(m);
        return false;
    }
    int base = (!list.empty() && lo == 0) ? 0 : 1; // a vertex 0 means 0..n-1
    for (size_t i = 0; i < list.size(); i++) {
        const Edge& e = list[i];
        if (e.from < base || e.from >= base + n || e.to < base || e.to >= base + n) {
            err = "line " + std::to_string(lineOf[i]) + ": vertices must be between " +
                  std::to_string(base) + " and " + std::to_string(base + n - 1);
            return false;
        }
    }
    g.reset(static_cast<int>(n), base);
    for (const Edge& e : list) g.add(e.from - base, e.to - base, e.cost);
    return true;
}

static bool saveGraph(const std::string& file, const Graph& g) {
    std::ofstream out(file.c_str());
    if (!out) return false;
    if (!g.note.empty()) out << "# " << g.note << '\n';
    out << g.n << ' ' << g.m() << '\n';
    for (const Edge& e : g.edges)
        out << g.label(e.from) << ' ' << g.label(e.to) << ' ' << e.cost << '\n';
    return static_cast<bool>(out);
}

// ---------------------------------------------------------------- generator

// xorshift32 (Marsaglia, 2003): the platform uses the same generator, so the same seed
// gives the same graph in the program and on the platform
struct Rng {
    uint32_t x;
    explicit Rng(uint32_t seed) : x(seed ? seed : 2463534242u) {}
    uint32_t next() {
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        return x;
    }
    int range(int a, int b) { // an integer in [a, b]
        uint32_t span = static_cast<uint32_t>(static_cast<long long>(b) - a + 1);
        return static_cast<int>(a + static_cast<long long>(next() % span));
    }
};

// Random graph with n vertices and m arcs, without loops or parallel arcs, with every
// vertex reachable from vertex 1: the first n - 1 arcs form a tree rooted at 1, then
// the arc list is shuffled. Costs are in [a, b]. If a < 0 and noNegCycle is true, the
// cost of the arc u -> v is c + p[u] - p[v], with c in [0, b] and p[v] in [0, -a]: on a
// cycle the p terms cancel, so no cycle has a negative cost.
static bool generate(Graph& g, int n, long long m, int a, int b, uint32_t seed,
                     bool noNegCycle, std::string& err) {
    long long most = static_cast<long long>(n) * (n - 1) / 2;
    if (n < 1 || n > 1000000) {
        err = "n must be between 1 and 10^6";
        return false;
    }
    if (m < n - 1 || m > most || m > 5000000) {
        err = "m must be between n - 1 = " + std::to_string(n - 1) + " and " +
              std::to_string(std::min(most, 5000000LL));
        return false;
    }
    if (a > b || a < -1000000 || b > 1000000) {
        err = "the costs must satisfy -10^6 <= a <= b <= 10^6";
        return false;
    }
    bool shift = noNegCycle && a < 0;
    if (shift && b < 0) {
        err = "without negative cycles, b must be >= 0";
        return false;
    }
    Rng rng(seed);
    std::vector<int> perm(n - 1), p(n, 0);
    for (int i = 0; i < n - 1; i++) perm[i] = i + 1;
    for (int i = n - 2; i > 0; i--) std::swap(perm[i], perm[rng.range(0, i)]);
    if (shift)
        for (int v = 0; v < n; v++) p[v] = rng.range(0, -a);
    bool small = static_cast<unsigned long long>(n) * n <= 50000000ULL;
    std::vector<bool> mat(small ? static_cast<size_t>(n) * n : 0, false);
    std::unordered_set<unsigned long long> seen;
    auto key = [n](int u, int v) { return static_cast<unsigned long long>(u) * n + v; };
    auto taken = [&](int u, int v) {
        return small ? static_cast<bool>(mat[key(u, v)]) : seen.count(key(u, v)) > 0;
    };
    std::vector<Edge> list;
    list.reserve(static_cast<size_t>(m));
    auto take = [&](int u, int v) {
        if (small) mat[key(u, v)] = true; else seen.insert(key(u, v));
        int c = shift ? rng.range(0, b) + p[u] - p[v] : rng.range(a, b);
        list.push_back(Edge{u, v, c});
    };
    for (int i = 0; i < n - 1; i++) { // the tree: perm[i] gets an earlier parent
        int r = rng.range(0, i);
        take(r == 0 ? 0 : perm[r - 1], perm[i]);
    }
    while (static_cast<long long>(list.size()) < m) {
        int u = rng.range(0, n - 1), v = rng.range(0, n - 1);
        if (u != v && !taken(u, v)) take(u, v);
    }
    for (int i = static_cast<int>(list.size()) - 1; i > 0; i--)
        std::swap(list[i], list[rng.range(0, i)]);
    g.reset(n, 1);
    for (const Edge& e : list) g.add(e.from, e.to, e.cost);
    std::ostringstream os;
    os << "random graph: n = " << n << ", m = " << m << ", costs in [" << a << ", " << b
       << "], seed " << seed << (shift ? ", no negative cycles" : "");
    g.note = os.str();
    return true;
}

// ---------------------------------------------------------------- output and timing

static std::string distStr(Dist d) {
    if (d == INF) return "inf";
    if (d == NEG_INF) return "-inf";
    return std::to_string(d);
}

// the time with 3 significant digits: "0.0123 ms", "1.23 ms", "123 ms"
static std::string timeStr(double ms) {
    std::ostringstream os;
    int digits = ms >= 100 ? 0 : ms >= 10 ? 1 : ms >= 1 ? 2 : ms >= 0.1 ? 3 : 4;
    os << std::fixed << std::setprecision(digits) << ms << " ms";
    return os.str();
}

static std::string pathStr(const Graph& g, const std::vector<int>& p) {
    std::string s;
    for (size_t i = 0; i < p.size(); i++)
        s += (i ? " -> " : "") + std::to_string(g.label(p[i]));
    return s;
}

static void printGraph(const Graph& g) {
    int neg = 0, cmin = INT_MAX, cmax = INT_MIN;
    for (const Edge& e : g.edges) {
        neg += e.cost < 0;
        cmin = std::min(cmin, e.cost);
        cmax = std::max(cmax, e.cost);
    }
    std::cout << " " << plural(g.n, "vertex", "vertices") << ", "
              << plural(g.m(), "arc", "arcs");
    if (g.m() > 0) std::cout << ", costs in [" << cmin << ", " << cmax << "]";
    std::cout << ", " << plural(neg, "negative-cost arc", "negative-cost arcs") << '\n';
    int shown = g.m() <= 200 ? g.n : std::min(g.n, 20);
    for (int u = 0; u < shown; u++) {
        std::cout << " " << std::setw(4) << g.label(u) << ":";
        if (g.adj[u].empty()) std::cout << " (no arcs)";
        for (const Arc& a : g.adj[u])
            std::cout << ' ' << g.label(a.to) << " (" << a.cost << ")";
        std::cout << '\n';
    }
    if (shown < g.n)
        std::cout << " ... (" << plural(g.n - shown, "vertex", "vertices")
                  << " not shown)\n";
}

// the cost of the cycle; between two neighbours the cheapest arc is taken
static Dist cycleCost(const Graph& g, const std::vector<int>& c) {
    Dist sum = 0;
    for (size_t i = 0; i + 1 < c.size(); i++) {
        Dist best = INF;
        for (const Arc& a : g.adj[c[i]])
            if (a.to == c[i + 1]) best = std::min(best, static_cast<Dist>(a.cost));
        sum += best;
    }
    return sum;
}

static void printResult(const Graph& g, const Result& r, int s) {
    int reach = 0, far = -1;
    for (int v = 0; v < g.n; v++) {
        if (r.dist[v] == INF) continue;
        reach++;
        if (r.dist[v] != NEG_INF && (far < 0 || r.dist[v] > r.dist[far])) far = v;
    }
    std::cout << " vertices reachable from " << g.label(s) << ": " << reach << " of "
              << g.n;
    if (far >= 0)
        std::cout << "; the farthest: " << g.label(far) << ", at distance "
                  << r.dist[far];
    std::cout << '\n';
    if (r.negCycle)
        std::cout << " negative cycle: " << pathStr(g, r.cycle) << " (cost "
                  << cycleCost(g, r.cycle) << "); the vertices reachable from it have "
                  << "distance -inf\n";
    int shown = g.n <= 30 ? g.n : 10;
    std::cout << " vertex distance predecessor path\n";
    for (int v = 0; v < shown; v++) {
        bool none = r.pred[v] < 0 || r.dist[v] == NEG_INF;
        std::cout << " " << std::setw(6) << g.label(v) << ' ' << std::setw(10)
                  << distStr(r.dist[v]) << ' ' << std::setw(12)
                  << (none ? std::string("-") : std::to_string(g.label(r.pred[v])))
                  << " " << pathStr(g, pathTo(r, v)) << '\n';
    }
    if (shown < g.n)
        std::cout << " ... (" << plural(g.n - shown, "vertex", "vertices")
                  << " not shown; option 8 prints the path to a vertex)\n";
}

typedef Result (*Algorithm)(const Graph&, int);
const char* const NAMES = {"Dijkstra, no heap", "Dijkstra, min-heap", "Bellman-Ford"};
const Algorithm ALGS = {dijkstraArray, dijkstraHeap, bellmanFord};

// the mean time of one run, in milliseconds; runs repeat until 0.2 s have passed
static double timeMs(Algorithm alg, const Graph& g, int s, Result* keep) {
    typedef std::chrono::steady_clock Clock;
    int runs = 0;
    double total = 0;
    do {
        Clock::time_point t0 = Clock::now();
        Result r = alg(g, s);
        total += std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
        if (runs++ == 0 && keep) *keep = r;
    } while (total < 200 && runs < 1000);
    return total / runs;
}

static int differences(const Result& a, const Result& b) {
    int k = 0;
    for (size_t v = 0; v < a.dist.size(); v++) k += a.dist[v] != b.dist[v];
    return k;
}

static std::string same(int diff) {
    if (diff == 0) return "the same distances";
    return plural(diff, "different distance", "different distances");
}

static void compareAll(const Graph& g, int s) {
    Result r = {Result(0), Result(0), Result(0)};
    double t;
    for (int i = 0; i < 3; i++) t[i] = timeMs(ALGS[i], g, s, &r[i]);
    std::cout << " algorithm time arcs examined"
              << " distances decreased\n";
    for (int i = 0; i < 3; i++) {
        std::cout << " " << std::left << std::setw(21) << NAMES[i] << std::right
                  << std::setw(11) << timeStr(t[i]) << std::setw(16) << r[i].relax
                  << std::setw(20) << r[i].updates;
        if (i == 2) std::cout << " (" << plural(r[i].passes, "pass", "passes") << ")";
        std::cout << '\n';
    }
    std::cout << " Dijkstra without and with a heap: "
              << same(differences(r, r)) << '\n';
    if (r.negCycle) {
        std::cout << " Bellman-Ford found a negative cycle: Dijkstra's distances mean "
                  << "nothing.\n";
        return;
    }
    int d = differences(r, r);
    std::cout << " Dijkstra and Bellman-Ford: " << same(d);
    if (d) std::cout << " (Dijkstra is wrong because of the negative-cost arcs)";
    std::cout << '\n';
}

// The doubling experiment: random graphs with n0, 2 n0, 4 n0, ... vertices and costs in
//. If T(2n) / T(n) approaches 2^b, the time grows roughly as n^b.
static void doubling(int n0, int count, int degree) {
    std::ofstream csv("doubling.csv");
    csv << "n,m,dijkstra_array_ms,dijkstra_heap_ms,bellman_ford_ms,passes\n";
    std::cout << " n m Dijkstra, no heap Dijkstra, heap "
              << " Bellman-Ford\n";
    double prev = {0, 0, 0};
    long long n = n0;
    for (int i = 0; i < count; i++, n *= 2) {
        long long m = std::max(n - 1, degree > 0 ? n * degree : n * (n - 1) / 4);
        Graph g;
        std::string err = "n too large";
        uint32_t seed = static_cast<uint32_t>(1 + i);
        bool ok = n <= 1000000;
        if (!ok || !generate(g, static_cast<int>(n), m, 1, 100, seed, false, err)) {
            std::cout << " n = " << n << ": stopped (" << err << ")\n";
            break;
        }
        std::cout << std::setw(9) << n << std::setw(11) << m << std::flush;
        Result bf(0);
        double t;
        for (int k = 0; k < 3; k++) {
            t[k] = timeMs(ALGS[k], g, 0, k == 2 ? &bf : nullptr);
            std::ostringstream cell;
            cell << timeStr(t[k]);
            if (prev[k] > 0) cell << " x" << std::fixed << std::setprecision(1)
                                  << t[k] / prev[k];
            std::cout << std::setw(18) << cell.str() << std::flush;
            prev[k] = t[k];
        }
        std::cout << " (" << plural(bf.passes, "pass", "passes") << ")\n";
        csv << n << ',' << m << ',' << t << ',' << t << ',' << t << ','
            << bf.passes << '\n';
    }
    std::cout << " xk = T(2n) / T(n). The results are also in doubling.csv.\n";
}

// ---------------------------------------------------------------- menu

enum Ask { OK, BAD, END }; // valid answer, invalid answer, end of input

static Ask askLine(const std::string& prompt, std::string& line) {
    std::cout << prompt;
    if (!std::getline(std::cin, line)) return END;
    size_t a = line.find_first_not_of(" \t\r"), b = line.find_last_not_of(" \t\r");
    line = a == std::string::npos ? "" : line.substr(a, b - a + 1);
    return OK;
}

// exactly one integer in [lo, hi]; *dflt is used for an empty line
static Ask askNum(const std::string& prompt, long long lo, long long hi, long long& x,
                  const long long* dflt = nullptr) {
    std::string line;
    if (askLine(prompt, line) == END) return END;
    if (line.empty() && dflt) { x = *dflt; return OK; }
    std::vector<long long> v;
    if (!parseNums(line, v) || v.size() != 1 || v < lo || v > hi) {
        std::cout << "Invalid value: enter one integer between " << lo << " and " << hi
                  << ".\n";
        return BAD;
    }
    x = v;
    return OK;
}

static Ask askYes(const std::string& prompt, bool& yes) {
    std::string line;
    if (askLine(prompt + " (y/n): ", line) == END) return END;
    if (line == "y" || line == "yes") { yes = true; return OK; }
    if (line == "n" || line == "no") { yes = false; return OK; }
    std::cout << "Answer y or n.\n";
    return BAD;
}

struct Session {
    Graph g;
    Result last = Result(0); // the result of the last run (options 5-7)
    int lastAlg = -1, lastSrc = 0;
};

static Ask askVertex(const Graph& g, const std::string& prompt, int& v) {
    long long x = 0;
    Ask a = askNum(prompt, g.base, g.base + g.n - 1, x);
    v = static_cast<int>(x) - g.base;
    return a;
}

static Ask optLoadSave(Session& S, bool load) {
    std::string file, err;
    if (askLine("File name: ", file) == END) return END;
    if (file.empty()) { std::cout << "Empty name.\n"; return BAD; }
    if (!load) {
        std::cout << (saveGraph(file, S.g) ? "Saved the graph to " + file + ".\n"
                                          : "The file cannot be written.\n");
        return OK;
    }
    Graph h;
    if (!loadGraph(file, h, err)) {
        std::cout << file << ": " << err << ".\n";
        return BAD;
    }
    S.g = h;
    S.lastAlg = -1;
    std::cout << "Loaded the graph.\n";
    printGraph(S.g);
    return OK;
}

static Ask optGenerate(Session& S) {
    long long n = 0, m = 0, a = 0, b = 0, seed = 0, one = 1;
    bool avoid = false;
    Ask r = askNum("Number of vertices n: ", 1, 1000000, n);
    if (r == OK) r = askNum("Number of arcs m: ", 0, 5000000, m);
    if (r == OK) r = askNum("Minimum cost a: ", -1000000, 1000000, a);
    if (r == OK) r = askNum("Maximum cost b: ", -1000000, 1000000, b);
    if (r == OK) r = askNum("Seed (Enter = 1): ", 0, UINT32_MAX, seed, &one);
    if (r == OK && a < 0) r = askYes("No negative cycles?", avoid);
    if (r != OK) return r;
    Graph h;
    std::string err;
    if (!generate(h, static_cast<int>(n), m, static_cast<int>(a), static_cast<int>(b),
                  static_cast<uint32_t>(seed), avoid, err)) {
        std::cout << "Invalid parameters: " << err << ".\n";
        return BAD;
    }
    S.g = h;
    S.lastAlg = -1;
    std::cout << "Generated a " << S.g.note << ".\n";
    return OK;
}

static Ask optRun(Session& S, int k) {
    int s = 0;
    Ask r = askVertex(S.g, "Source vertex: ", s);
    if (r != OK) return r;
    if (k < 2 && S.g.hasNegative())
        std::cout << "Warning: the graph has negative-cost arcs, so Dijkstra may give "
                  << "wrong distances. Use Bellman-Ford.\n";
    double t = timeMs(ALGS[k], S.g, s, &S.last);
    S.lastAlg = k;
    S.lastSrc = s;
    std::cout << NAMES[k] << " from vertex " << S.g.label(s) << ": " << timeStr(t)
              << ", " << plural(S.last.relax, "arc examined", "arcs examined") << ", "
              << plural(S.last.updates, "distance decreased", "distances decreased");
    if (k == 2) std::cout << ", " << plural(S.last.passes, "pass", "passes");
    std::cout << ".\n";
    printResult(S.g, S.last, s);
    return OK;
}

static Ask optPath(Session& S) {
    if (S.lastAlg < 0) {
        std::cout << "Run one of the algorithms first (options 5-7).\n";
        return BAD;
    }
    int t = 0;
    Ask r = askVertex(S.g, "Target vertex: ", t);
    if (r != OK) return r;
    std::cout << " " << NAMES[S.lastAlg] << ", from " << S.g.label(S.lastSrc) << " to "
              << S.g.label(t) << ": ";
    std::vector<int> p = pathTo(S.last, t);
    if (S.last.dist[t] == INF)
        std::cout << "there is no path.\n";
    else if (S.last.dist[t] == NEG_INF)
        std::cout << "distance -inf: a path can go around the negative cycle any "
                  << "number of times.\n";
    else
        std::cout << "distance " << S.last.dist[t] << ", path " << pathStr(S.g, p)
                  << " (" << plural(static_cast<long long>(p.size()) - 1, "arc", "arcs")
                  << ").\n";
    return OK;
}

static Ask optDoubling() {
    long long n0 = 0, count = 0, degree = 0, d0 = 1000, c0 = 5, g0 = 5;
    Ask r = askNum("Initial number of vertices n0 (Enter = 1000): ", 2, 1000000, n0, &d0);
    if (r == OK) r = askNum("Number of graphs (Enter = 5): ", 1, 12, count, &c0);
    if (r == OK)
        r = askNum("Mean degree m / n, 0 for a dense graph (Enter = 5): ", 0, 1000,
                   degree, &g0);
    if (r != OK) return r;
    doubling(static_cast<int>(n0), static_cast<int>(count), static_cast<int>(degree));
    return OK;
}

int main(int argc, char** argv) {
    Session S;
    std::string err;
    std::cout << "Single-source shortest paths - DSA lab\n";
    if (argc > 1) {
        if (loadGraph(argv, S.g, err))
            std::cout << "Loaded " << argv << ".\n";
        else
            std::cout << argv << ": " << err << ".\n";
    }
    for (;;) {
        std::cout << "\n=== graph: " << plural(S.g.n, "vertex", "vertices") << ", "
                  << plural(S.g.m(), "arc", "arcs") << " ===\n"
                  << " 1. Load a graph from a file\n"
                  << " 2. Generate a random graph\n"
                  << " 3. Save the graph to a file\n"
                  << " 4. Print the graph\n"
                  << " 5. Dijkstra without a heap, O(n^2 + m)\n"
                  << " 6. Dijkstra with a min-heap, O((n + m) log n)\n"
                  << " 7. Bellman-Ford, O(n m)\n"
                  << " 8. Shortest path to a vertex, after options 5-7\n"
                  << " 9. Compare the three algorithms on the current graph\n"
                  << "10. The doubling experiment, on random graphs\n"
                  << " 0. Exit\n";
        long long o = 0;
        Ask r = askNum("Option: ", 0, 10, o);
        if (r == END) break;
        if (r == BAD) continue;
        if (o == 0) break;
        if (o >= 3 && o <= 9 && S.g.n == 0) {
            std::cout << "There is no graph yet: load one (1) or generate one (2).\n";
            continue;
        }
        int s = 0;
        if (o == 1 || o == 3) {
            r = optLoadSave(S, o == 1);
        } else if (o == 2) {
            r = optGenerate(S);
        } else if (o == 4) {
            printGraph(S.g);
        } else if (o <= 7) {
            r = optRun(S, static_cast<int>(o - 5));
        } else if (o == 8) {
            r = optPath(S);
        } else if (o == 9) {
            r = askVertex(S.g, "Source vertex: ", s);
            if (r == OK) compareAll(S.g, s);
        } else {
            r = optDoubling();
        }
        if (r == END) break;
    }
    std::cout << "Goodbye!\n";
    return 0;
}