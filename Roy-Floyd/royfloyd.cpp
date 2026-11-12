// All-pairs shortest paths - the Roy-Floyd (Floyd-Warshall) algorithm - DSA lab
//
// For a directed graph with arc costs, the program computes the distance matrix
// D, in which D[i][j] is the cost of a shortest path from i to j, and the
// matrix P, in which P[i][j] is the vertex before j on a shortest path from i
// to j:
// - Roy-Floyd (Roy-Floyd-Warshall): n iterations, one for every vertex k; in
// iteration k, paths may also go through k, O(n^3);
// - Warshall: the transitive closure, the same scheme with Boolean values;
// - n runs of Dijkstra, without a heap and with a min-heap, for comparison.
// Roy-Floyd accepts negative costs; a negative-cost cycle shows up on the
// diagonal, as D[v][v] < 0. Dijkstra gives correct distances only if no arc has
// a negative cost.
//
// A graph file: n and m on the first line, then m lines "u v c", the arc u -> v
// of cost c. The arcs may be followed by lines "x y", the pairs whose shortest
// path is asked for (option 7). Vertices are numbered from 1 to n; the
// numbering 0..n-1 is accepted as well. Lines that start with # are comments.
//
// Build: g++ -std=c++11 -O2 -Wall -Wextra -o roy_floyd roy_floyd_en.cpp
// Run: ./roy_floyd or ./roy_floyd graf_exemplu.txt

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

typedef long long Dist;
const Dist INF = LLONG_MAX / 4;  // there is no path (yet)
const Dist NEG_INF = -INF;       // a path that can go through a negative cycle
const int MAX_COST = 1000000000; // |c| <= 10^9, so that sums fit in long long
const int MAX_N = 2000;          // D has n^2 elements: 32 MB for n = 2000

struct Arc {
  int to, cost;
}; // arc in the adjacency list of a vertex
struct Edge {
  int from, to, cost;
}; // arc in the list of all the arcs

struct Graph {
  int n = 0;               // the vertices are 0..n-1 in memory
  int base = 1;            // the numbering of the file: from 1 or from 0
  std::string note;        // comment written when saving (from the generator)
  std::vector<Edge> edges; // the arcs, in the order of the file
  std::vector<std::vector<Arc>> adj;        // the adjacency lists (Dijkstra)
  std::vector<std::pair<int, int>> queries; // the pairs "x y" after the arcs

  void reset(int vertices, int numbering) {
    n = vertices;
    base = numbering;
    note.clear();
    edges.clear();
    queries.clear();
    adj.assign(n, std::vector<Arc>());
  }

  void add(int u, int v, int c) {
    edges.push_back(Edge{u, v, c});
    adj[u].push_back(Arc{v, c});
  }

  int m() const { return static_cast<int>(edges.size()); }
  int label(int v) const { return v + base; }

  bool hasNegative() const {
    for (const Edge &e : edges)
      if (e.cost < 0)
        return true;
    return false;
  }
};

// an n x n matrix, stored row after row in a single vector: a[i][j] is the
// element in row i and column j, and a whole row sits in memory in one block
template <typename T> struct Matrix {
  int n = 0;
  std::vector<T> a;
  Matrix() {}
  Matrix(int size, T value)
      : n(size), a(static_cast<size_t>(size) * size, value) {}
  T *operator[](int i) { return &a[static_cast<size_t>(i) * n]; }
  const T *operator[](int i) const { return &a[static_cast<size_t>(i) * n]; }
};

// D[i][j]: the distance from i to j (INF: no path; NEG_INF: -inf); P[i][j]: the
// vertex before j on the shortest path from i, or -1
struct Apsp {
  Matrix<Dist> d;
  Matrix<int> p;
  long long compared = 0; // Roy-Floyd: comparisons; Dijkstra: arcs examined
  long long updates = 0;  // values decreased
  bool negCycle = false;  // there is a negative-cost cycle
  std::vector<int>
      onCycle; // the vertices v with D[v][v] < 0, on negative cycles
};

// ---------------------------------------------------------------- Roy-Floyd

// D^(0) and P^(0): 0 on the diagonal, the cost of the arc u -> v (the smallest
// one, if there are several arcs u -> v) and INF elsewhere; P[u][v] = u for
// every arc. A loop u -> u of negative cost is already a negative cycle:
// D[u][u] = c < 0.
static void initApsp(const Graph &g, Apsp &r) {
  r.d = Matrix<Dist>(g.n, INF);
  r.p = Matrix<int>(g.n, -1);
  for (int v = 0; v < g.n; v++)
    r.d[v][v] = 0;
  for (const Edge &e : g.edges)
    if (e.cost < r.d[e.from][e.to]) {
      r.d[e.from][e.to] = e.cost;
      r.p[e.from][e.to] = e.from;
    }
}

// with a negative cycle: D[i][j] = -inf if a vertex v with D[v][v] < 0 is
// reachable from i and j is reachable from v; such a path can go around the
// cycle any number of times
static void markNegInf(Apsp &r) {
  int n = r.d.n;
  for (int v = 0; v < n; v++)
    if (r.d[v][v] < 0)
      r.onCycle.push_back(v);
  if (r.onCycle.empty())
    return;
  r.negCycle = true;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      for (int v : r.onCycle)
        if (r.d[i][v] != INF && r.d[v][j] != INF) {
          r.d[i][j] = NEG_INF;
          r.p[i][j] = -1;
          break;
        }
}

// Roy-Floyd: after iteration k, D[i][j] is the cost of a shortest path from i
// to j whose intermediate vertices are all among the first k vertices. A path
// that also goes through vertex k is made of two paths already computed, i -> k
// and k -> j. Row k and column k do not change in iteration k, so the matrix
// can be updated in place.
Apsp royFloyd(const Graph &g) {
  Apsp r;
  initApsp(g, r);
  int n = g.n;
  // local counters: a counter inside r, increased in the inner loop, would
  // force the compiler to write it to memory at every step (di and dk might
  // point to it)
  long long compared = 0, updates = 0;
  for (int k = 0; k < n; k++) {
    const Dist *dk = r.d[k];
    const int *pk = r.p[k];
    for (int i = 0; i < n; i++) {
      Dist dik = r.d[i][k];
      if (dik == INF)
        continue; // no path i -> k: row i does not change
      Dist *di = r.d[i];
      int *pi = r.p[i];
      for (int j = 0; j < n; j++) {
        if (dk[j] == INF)
          continue; // no path k -> j
        compared++;
        if (dik + dk[j] < di[j]) {
          // max: with a negative cycle, the values can drop without bound
          di[j] = std::max(dik + dk[j], NEG_INF);
          pi[j] = pk[j]; // the path k -> j ends the same way
          updates++;
        }
      }
    }
  }
  r.compared = compared;
  r.updates = updates;
  markNegInf(r);
  return r;
}

// Warshall: the transitive closure. T[i][j] = 1 if j is reachable from i (every
// vertex is reachable from itself). The same scheme as Roy-Floyd: j becomes
// reachable from i through k if k is reachable from i and j is reachable from
// k.
Matrix<char> warshall(const Graph &g, long long &updates) {
  int n = g.n;
  Matrix<char> t(n, 0);
  for (int v = 0; v < n; v++)
    t[v][v] = 1;
  for (const Edge &e : g.edges)
    t[e.from][e.to] = 1;
  updates = 0;
  for (int k = 0; k < n; k++)
    for (int i = 0; i < n; i++) {
      if (!t[i][k])
        continue;
      for (int j = 0; j < n; j++)
        if (t[k][j] && !t[i][j]) {
          t[i][j] = 1;
          updates++;
        }
    }
  return t;
}

// the shortest path from s to t, rebuilt from P: the vertex before t is
// P[s][t], the one before P[s][t] is P[s][P[s][t]] and so on, back to s; empty
// if there is no finite path
static std::vector<int> pathP(const Apsp &r, int s, int t) {
  std::vector<int> p;
  if (r.d[s][t] == INF || r.d[s][t] == NEG_INF)
    return p;
  for (int v = t; v != s; v = r.p[s][v]) {
    if (v < 0 || static_cast<int>(p.size()) > r.d.n)
      return std::vector<int>();
    p.push_back(v);
  }
  p.push_back(s);
  std::reverse(p.begin(), p.end());
  return p;
}

// The method of the ADS implementation on GitHub, without the matrix P: if
// D[i][j] = D[i][k] + D[k][j] for a vertex k other than i and j, a shortest
// path goes through k and is rebuilt in pieces: i -> k and k -> j. Otherwise,
// the shortest path is the arc i -> j itself. Every split searches for k in
// O(n). With a cycle of cost 0, the splits can repeat forever, so the depth and
// the length of the path are limited to n.
static bool splitPath(const Apsp &r, int i, int j, int depth,
                      std::vector<int> &out) {
  int n = r.d.n;
  if (depth > n || static_cast<int>(out.size()) > n)
    return false;
  for (int k = 0; k < n; k++)
    if (k != i && k != j && r.d[i][k] != INF && r.d[k][j] != INF &&
        r.d[i][k] + r.d[k][j] == r.d[i][j])
      return splitPath(r, i, k, depth + 1, out) &&
             splitPath(r, k, j, depth + 1, out);
  out.push_back(j);
  return true;
}

static std::vector<int> pathSplit(const Apsp &r, int s, int t) {
  std::vector<int> p;
  if (r.d[s][t] == INF || r.d[s][t] == NEG_INF)
    return p;
  p.push_back(s);
  if (s != t && !splitPath(r, s, t, 0, p))
    return std::vector<int>();
  return p;
}

// ------------------------------------------------------- Dijkstra from every
// vertex

// Min-heap of vertices ordered by (d[v], v), as in the Dijkstra lab. pos[v] is
// the position of v in the heap or -1, so that decrease finds the vertex in
// O(1).
class MinHeap {
public:
  MinHeap(int n, const std::vector<Dist> &d) : d_(d), pos_(n, -1) {
    h_.reserve(n);
  }
  bool empty() const { return h_.empty(); }
  bool contains(int v) const { return pos_[v] >= 0; }

  void push(int v) {
    h_.push_back(v);
    up(static_cast<int>(h_.size()) - 1);
  }

  void decrease(int v) { up(pos_[v]); }

  int pop() {
    int top = h_.front(), last = h_.back();
    h_.pop_back();
    pos_[top] = -1;
    if (!h_.empty()) {
      h_[0] = last;
      pos_[last] = 0;
      down(0);
    }
    return top;
  }

private:
  const std::vector<Dist> &d_;
  std::vector<int> h_, pos_;

  bool less(int a, int b) const {
    return d_[a] < d_[b] || (d_[a] == d_[b] && a < b);
  }
  void place(int v, int i) {
    h_[i] = v;
    pos_[v] = i;
  }

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
      if (c >= size)
        break;
      if (c + 1 < size && less(h_[c + 1], h_[c]))
        c++;
      if (!less(h_[c], v))
        break;
      place(h_[c], i);
      i = c;
    }
    place(v, i);
  }
};

// n runs of Dijkstra, one from every vertex s; row s of D is the array d of the
// run from s, and row s of P is the array pred. Without a heap: O(n (n^2 + m))
// = O(n^3); with a min-heap: O(n (n + m) log n). Correct only without
// negative-cost arcs.
Apsp dijkstraAll(const Graph &g, bool useHeap) {
  Apsp r;
  int n = g.n;
  r.d = Matrix<Dist>(n, INF);
  r.p = Matrix<int>(n, -1);
  std::vector<Dist> dist(n);
  std::vector<int> pred(n);
  std::vector<char> done(n);

  for (int s = 0; s < n; s++) {
    std::fill(dist.begin(), dist.end(), INF);
    std::fill(pred.begin(), pred.end(), -1);
    std::fill(done.begin(), done.end(), 0);
    dist[s] = 0;
    MinHeap h(n, dist);
    if (useHeap)
      h.push(s);
    for (;;) {
      int u = -1;
      if (useHeap) {
        if (h.empty())
          break;
        u = h.pop();
      } else {
        for (int v = 0; v < n; v++)
          if (!done[v] && dist[v] < INF && (u < 0 || dist[v] < dist[u]))
            u = v;
        if (u < 0)
          break;
      }
      done[u] = 1;
      for (const Arc &a : g.adj[u]) {
        r.compared++;
        if (!done[a.to] && dist[u] + a.cost < dist[a.to]) {
          dist[a.to] = dist[u] + a.cost;
          pred[a.to] = u;
          r.updates++;
          if (!useHeap)
            continue;
          if (h.contains(a.to))
            h.decrease(a.to);
          else
            h.push(a.to);
        }
      }
    }
    std::copy(dist.begin(), dist.end(), r.d[s]);
    std::copy(pred.begin(), pred.end(), r.p[s]);
  }
  return r;
}

// ---------------------------------------------------------------- files

// "1 arc", "5 arcs"
static std::string plural(long long n, const char *one, const char *many) {
  if (n == 1)
    return std::string("1 ") + one;
  return std::to_string(n) + " " + many;
}

static bool isSep(char c) {
  return c == ' ' || c == '\t' || c == ',' || c == ';' || c == '\r';
}

// reads the integers of a line; separators: spaces, commas, semicolons
static bool parseNums(const std::string &line, std::vector<long long> &out) {
  out.clear();
  size_t i = 0;
  while (i < line.size()) {
    if (isSep(line[i])) {
      i++;
      continue;
    }
    size_t j = i;
    while (j < line.size() && !isSep(line[j]))
      j++;
    std::string tok = line.substr(i, j - i);
    errno = 0;
    char *end = nullptr;
    long long v = std::strtoll(tok.c_str(), &end, 10);
    if (*end != '\0' || errno == ERANGE)
      return false;
    out.push_back(v);
    i = j;
  }
  return true;
}

static bool blankOrComment(const std::string &line) {
  for (char c : line) {
    if (c == '#')
      return true;
    if (!isSep(c))
      return false;
  }
  return true;
}

// reads the graph from a file; on error returns false and explains in err
static bool loadGraph(const std::string &file, Graph &g, std::string &err) {
  std::ifstream in(file.c_str());
  if (!in) {
    err = "the file cannot be opened";
    return false;
  }
  std::string line;
  std::vector<long long> v;
  long long n = -1, m = -1;
  std::vector<Edge> list;
  std::vector<std::pair<long long, long long>> pairs;
  std::vector<int> lineOf, pairLine;
  int no = 0;
  long long lo = LLONG_MAX;
  while (std::getline(in, line)) {
    no++;
    if (blankOrComment(line))
      continue;
    std::string where = "line " + std::to_string(no) + ": ";
    if (!parseNums(line, v)) {
      err = where + "integers were expected";
      return false;
    }
    if (n < 0) {
      if (v.size() != 2) {
        err = where + "n and m were expected";
        return false;
      }
      n = v[0];
      m = v[1];
      if (n < 1 || n > MAX_N) {
        err = where + "n must be between 1 and " + std::to_string(MAX_N);
        return false;
      }
      if (m < 0 || m > 5000000) {
        err = where + "m must be between 0 and 5*10^6";
        return false;
      }
      list.reserve(static_cast<size_t>(m));
      continue;
    }
    if (static_cast<long long>(list.size()) == m) { // after the arcs: pairs x y
      if (v.size() != 2) {
        err = where +
              (v.size() == 3
                   ? "the file has more than m = " + std::to_string(m) + " arcs"
                   : std::string("a pair x y was expected"));
        return false;
      }
      pairs.push_back(std::make_pair(v[0], v[1]));
      pairLine.push_back(no);
      lo = std::min(lo, std::min(v[0], v[1]));
      continue;
    }
    if (v.size() != 3) {
      err = where + "an arc u v c was expected";
      return false;
    }
    if (v[0] < INT_MIN || v[0] > INT_MAX || v[1] < INT_MIN || v[1] > INT_MAX) {
      err = where + "vertex out of range";
      return false;
    }
    if (v[2] < -MAX_COST || v[2] > MAX_COST) {
      err = where + "the cost must be between -10^9 and 10^9";
      return false;
    }
    Edge e = {static_cast<int>(v[0]), static_cast<int>(v[1]),
              static_cast<int>(v[2])};
    lo = std::min(lo, static_cast<long long>(std::min(e.from, e.to)));
    list.push_back(e);
    lineOf.push_back(no);
  }
  if (n < 0) {
    err = "the file does not hold n and m";
    return false;
  }
  if (static_cast<long long>(list.size()) != m) {
    long long k = static_cast<long long>(list.size());
    err = "the file has " + plural(k, "arc", "arcs") +
          ", not m = " + std::to_string(m);
    return false;
  }
  int base = lo == 0 ? 0 : 1; // a vertex 0 means the numbering 0..n-1
  std::string range = "vertices must be between " + std::to_string(base) +
                      " and " + std::to_string(base + n - 1);
  for (size_t i = 0; i < list.size(); i++) {
    const Edge &e = list[i];
    if (e.from < base || e.from >= base + n || e.to < base ||
        e.to >= base + n) {
      err = "line " + std::to_string(lineOf[i]) + ": " + range;
      return false;
    }
  }
  for (size_t i = 0; i < pairs.size(); i++)
    if (pairs[i].first < base || pairs[i].first >= base + n ||
        pairs[i].second < base || pairs[i].second >= base + n) {
      err = "line " + std::to_string(pairLine[i]) + ": " + range;
      return false;
    }
  g.reset(static_cast<int>(n), base);
  for (const Edge &e : list)
    g.add(e.from - base, e.to - base, e.cost);
  for (const auto &q : pairs)
    g.queries.push_back(std::make_pair(static_cast<int>(q.first) - base,
                                       static_cast<int>(q.second) - base));
  return true;
}

static bool saveGraph(const std::string &file, const Graph &g) {
  std::ofstream out(file.c_str());
  if (!out)
    return false;
  if (!g.note.empty())
    out << "# " << g.note << '\n';
  out << g.n << ' ' << g.m() << '\n';
  for (const Edge &e : g.edges)
    out << g.label(e.from) << ' ' << g.label(e.to) << ' ' << e.cost << '\n';
  for (const auto &q : g.queries)
    out << g.label(q.first) << ' ' << g.label(q.second) << '\n';
  return static_cast<bool>(out);
}

// ---------------------------------------------------------------- generator

// xorshift32 (Marsaglia, 2003): the platform uses the same generator, so the
// same seed gives the same graph in the program and on the platform
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

// Random graph with n vertices and m arcs, without loops or parallel arcs, with
// every vertex reachable from vertex 1: the first n - 1 arcs form a tree rooted
// at 1, then the arc list is shuffled. Costs are in [a, b]. If a < 0 and
// noNegCycle is true, the cost of the arc u -> v is c + p[u] - p[v], with c in
// [0, b] and p[v] in [0, -a]: on a cycle the p terms cancel, so no cycle has a
// negative cost. It is the same generator as in the Dijkstra and Bellman-Ford
// lab.
static bool generate(Graph &g, int n, long long m, int a, int b, uint32_t seed,
                     bool noNegCycle, std::string &err) {
  long long most = static_cast<long long>(n) * (n - 1) / 2;
  if (n < 1 || n > MAX_N) {
    err = "n must be between 1 and " + std::to_string(MAX_N);
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
  for (int i = 0; i < n - 1; i++)
    perm[i] = i + 1;
  for (int i = n - 2; i > 0; i--)
    std::swap(perm[i], perm[rng.range(0, i)]);
  if (shift)
    for (int v = 0; v < n; v++)
      p[v] = rng.range(0, -a);
  bool small = static_cast<unsigned long long>(n) * n <= 50000000ULL;
  std::vector<bool> mat(small ? static_cast<size_t>(n) * n : 0, false);
  std::unordered_set<unsigned long long> seen;
  auto key = [n](int u, int v) {
    return static_cast<unsigned long long>(u) * n + v;
  };
  auto taken = [&](int u, int v) {
    return small ? static_cast<bool>(mat[key(u, v)])
                 : seen.count(key(u, v)) > 0;
  };
  std::vector<Edge> list;
  list.reserve(static_cast<size_t>(m));
  auto take = [&](int u, int v) {
    if (small)
      mat[key(u, v)] = true;
    else
      seen.insert(key(u, v));
    int c = shift ? rng.range(0, b) + p[u] - p[v] : rng.range(a, b);
    list.push_back(Edge{u, v, c});
  };
  for (int i = 0; i < n - 1;
       i++) { // the tree: the parent of perm[i] comes before it
    int r = rng.range(0, i);
    take(r == 0 ? 0 : perm[r - 1], perm[i]);
  }
  while (static_cast<long long>(list.size()) < m) {
    int u = rng.range(0, n - 1), v = rng.range(0, n - 1);
    if (u != v && !taken(u, v))
      take(u, v);
  }
  for (int i = static_cast<int>(list.size()) - 1; i > 0; i--)
    std::swap(list[i], list[rng.range(0, i)]);
  g.reset(n, 1);
  for (const Edge &e : list)
    g.add(e.from, e.to, e.cost);
  std::ostringstream os;
  os << "random graph: n = " << n << ", m = " << m << ", costs in [" << a
     << ", " << b << "], seed " << seed
     << (shift ? ", no negative cycles" : "");
  g.note = os.str();
  return true;
}

// ---------------------------------------------------------------- output and
// timing

static std::string distStr(Dist d) {
  if (d == INF)
    return "-"; // no path, as in the original program
  if (d == NEG_INF)
    return "-inf";
  return std::to_string(d);
}

// the time with 3 significant digits: "0.0123 ms", "1.23 ms", "123 ms"
static std::string timeStr(double ms) {
  std::ostringstream os;
  int digits = ms >= 100 ? 0 : ms >= 10 ? 1 : ms >= 1 ? 2 : ms >= 0.1 ? 3 : 4;
  os << std::fixed << std::setprecision(digits) << ms << " ms";
  return os.str();
}

static std::string pathStr(const Graph &g, const std::vector<int> &p) {
  std::string s;
  for (size_t i = 0; i < p.size(); i++)
    s += (i ? " -> " : "") + std::to_string(g.label(p[i]));
  return s;
}

// an n x n matrix with the column header; mark[i][j] != 0 adds a star after the
// value
static void printMatrix(const Graph &g, const std::vector<std::string> &cell,
                        const std::vector<char> *mark) {
  int n = g.n;
  size_t w = std::max<size_t>(3, std::to_string(g.label(n - 1)).size());
  for (const std::string &s : cell)
    w = std::max(w, s.size());
  w += 2;
  std::cout << " " << std::setw(static_cast<int>(w)) << ' ';
  for (int j = 0; j < n; j++)
    std::cout << std::setw(static_cast<int>(w)) << g.label(j) << ' ';
  std::cout << '\n';
  for (int i = 0; i < n; i++) {
    std::cout << " " << std::setw(static_cast<int>(w)) << g.label(i);
    for (int j = 0; j < n; j++) {
      size_t k = static_cast<size_t>(i) * n + j;
      std::cout << std::setw(static_cast<int>(w)) << cell[k]
                << (mark && (*mark)[k] ? '*' : ' ');
    }
    std::cout << '\n';
  }
}

static void printDist(const Graph &g, const Matrix<Dist> &d,
                      const std::vector<char> *mark) {
  std::vector<std::string> cell;
  for (Dist x : d.a)
    cell.push_back(distStr(x));
  printMatrix(g, cell, mark);
}

static void printPred(const Graph &g, const Matrix<int> &p) {
  std::vector<std::string> cell;
  for (int x : p.a)
    cell.push_back(x < 0 ? "-" : std::to_string(g.label(x)));
  printMatrix(g, cell, nullptr);
}

static void printGraph(const Graph &g) {
  int neg = 0, cmin = INT_MAX, cmax = INT_MIN;
  for (const Edge &e : g.edges) {
    neg += e.cost < 0;
    cmin = std::min(cmin, e.cost);
    cmax = std::max(cmax, e.cost);
  }
  std::cout << " " << plural(g.n, "vertex", "vertices") << ", "
            << plural(g.m(), "arc", "arcs");
  if (g.m() > 0)
    std::cout << ", costs in [" << cmin << ", " << cmax << "]";
  std::cout << ", " << neg << " with a negative cost\n";
  int shown = g.m() <= 200 ? g.n : std::min(g.n, 20);
  for (int u = 0; u < shown; u++) {
    std::cout << " " << std::setw(4) << g.label(u) << ":";
    if (g.adj[u].empty())
      std::cout << " (no arcs)";
    for (const Arc &a : g.adj[u])
      std::cout << ' ' << g.label(a.to) << " (" << a.cost << ")";
    std::cout << '\n';
  }
  if (shown < g.n)
    std::cout << " ... (" << plural(g.n - shown, "vertex", "vertices")
              << " not shown)\n";
  if (g.n <= 20) {
    Apsp r;
    initApsp(g, r);
    std::cout << " The cost matrix, D^(0) (\"-\" = no arc):\n";
    printDist(g, r.d, nullptr);
  }
  if (!g.queries.empty())
    std::cout << " The file asks for "
              << plural(static_cast<long long>(g.queries.size()), "path",
                        "paths")
              << " (option 7).\n";
}

static std::string vertexList(const Graph &g, const std::vector<int> &vs) {
  std::string s;
  for (size_t i = 0; i < vs.size(); i++)
    s += (i ? ", " : "") + std::to_string(g.label(vs[i]));
  return s;
}

static void printResult(const Graph &g, const Apsp &r) {
  long long reach = 0, neg = 0;
  int a = -1, b = -1;
  for (int i = 0; i < g.n; i++)
    for (int j = 0; j < g.n; j++) {
      Dist x = r.d[i][j];
      if (i == j || x == INF)
        continue;
      reach++;
      if (x == NEG_INF) {
        neg++;
        continue;
      }
      if (a < 0 || x > r.d[a][b]) {
        a = i;
        b = j;
      }
    }
  long long all = static_cast<long long>(g.n) * (g.n - 1);
  std::cout << " pairs (i, j), i != j, with a path from i to j: " << reach
            << " of " << all;
  if (a >= 0)
    std::cout << "; the largest distance: D[" << g.label(a) << "]["
              << g.label(b) << "] = " << r.d[a][b];
  std::cout << '\n';
  if (r.negCycle)
    std::cout << " negative cycles: D[v][v] < 0 for v = "
              << vertexList(g, r.onCycle) << "; "
              << plural(neg, "pair has", "pairs have") << " distance -inf\n";
  if (g.n <= 20) {
    std::cout << " The distance matrix D (\"-\" = no path):\n";
    printDist(g, r.d, nullptr);
    std::cout << " The predecessor matrix P (P[i][j] = the vertex before j):\n";
    printPred(g, r.p);
  } else {
    std::cout << " (the matrices have " << g.n << " x " << g.n
              << " elements and are "
              << "not shown; option 7 shows the path between two vertices)\n";
  }
}

// the same loop as royFloyd, printing every value that decreases and the matrix
// D after every iteration k; the values changed in iteration k are marked with
// *
static void traceRoyFloyd(const Graph &g) {
  Apsp r;
  initApsp(g, r);
  int n = g.n;
  std::cout << " D^(0), the cost matrix:\n";
  printDist(g, r.d, nullptr);
  for (int k = 0; k < n; k++) {
    std::vector<char> mark(static_cast<size_t>(n) * n, 0);
    int changed = 0;
    std::cout << "\n k = " << g.label(k)
              << ": paths may also go through vertex " << g.label(k) << '\n';
    for (int i = 0; i < n; i++) {
      Dist dik = r.d[i][k];
      if (dik == INF)
        continue;
      for (int j = 0; j < n; j++) {
        Dist dkj = r.d[k][j];
        if (dkj == INF || dik + dkj >= r.d[i][j])
          continue;
        std::cout << " D[" << g.label(i) << "][" << g.label(j)
                  << "]: " << distStr(r.d[i][j]) << " -> " << dik + dkj
                  << " (D[" << g.label(i) << "][" << g.label(k) << "] + D["
                  << g.label(k) << "][" << g.label(j) << "] = " << dik << " + "
                  << dkj << ")\n";
        r.d[i][j] = std::max(dik + dkj, NEG_INF);
        r.p[i][j] = r.p[k][j];
        mark[static_cast<size_t>(i) * n + j] = 1;
        changed++;
      }
    }
    if (!changed)
      std::cout << " no value changes\n";
    std::cout << " D^(" << g.label(k) << "):\n";
    printDist(g, r.d, &mark);
  }
  for (int v = 0; v < n; v++)
    if (r.d[v][v] < 0) {
      std::cout << " D[" << g.label(v) << "][" << g.label(v) << "] < 0: vertex "
                << g.label(v) << " lies on a negative-cost cycle.\n";
    }
}

typedef Apsp (*Algorithm)(const Graph &);
static Apsp dijkstraArrayAll(const Graph &g) { return dijkstraAll(g, false); }
static Apsp dijkstraHeapAll(const Graph &g) { return dijkstraAll(g, true); }
const char *const NAMES[3] = {"Roy-Floyd", "n x Dijkstra, no heap",
                              "n x Dijkstra, min-heap"};
const Algorithm ALGS[3] = {royFloyd, dijkstraArrayAll, dijkstraHeapAll};

// the mean time of one run, in milliseconds; runs repeat until 0.2 s have
// passed
static double timeMs(Algorithm alg, const Graph &g, Apsp *keep) {
  typedef std::chrono::steady_clock Clock;
  int runs = 0;
  double total = 0;
  do {
    Clock::time_point t0 = Clock::now();
    Apsp r = alg(g);
    total +=
        std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
    if (runs++ == 0 && keep)
      *keep = r;
  } while (total < 200 && runs < 1000);
  return total / runs;
}

static long long differences(const Apsp &a, const Apsp &b) {
  long long k = 0;
  for (size_t i = 0; i < a.d.a.size(); i++)
    k += a.d.a[i] != b.d.a[i];
  return k;
}

static void compareAll(const Graph &g) {
  Apsp r[3];
  double t[3];
  for (int i = 0; i < 3; i++)
    t[i] = timeMs(ALGS[i], g, &r[i]);
  std::cout << " algorithm time comparisons"
            << " values decreased\n";
  for (int i = 0; i < 3; i++)
    std::cout << " " << std::left << std::setw(23) << NAMES[i] << std::right
              << std::setw(13) << timeStr(t[i]) << std::setw(16)
              << r[i].compared << std::setw(18) << r[i].updates << '\n';
  std::cout << " (for Dijkstra, the comparisons are the arcs examined)\n";
  if (r[0].negCycle) {
    std::cout << " The graph has a negative cycle: Dijkstra's distances mean "
              << "nothing.\n";
    return;
  }
  long long d1 = differences(r[0], r[1]), d2 = differences(r[0], r[2]);
  if (!d1 && !d2)
    std::cout << " All three algorithms give the same matrix D.\n";
  else
    std::cout << " Roy-Floyd and Dijkstra: "
              << plural(std::max(d1, d2), "distance differs",
                        "distances differ")
              << " (Dijkstra is wrong because of the negative-cost arcs)\n";
}

// The doubling experiment: random graphs with n0, 2 n0, 4 n0, ... vertices and
// costs in
//. If T(2n) / T(n) approaches 2^b, the time grows roughly as n^b.
static void doubling(int n0, int count, int degree) {
  std::ofstream csv("dublare_rf.csv");
  csv << "n,m,roy_floyd_ms,dijkstra_tablou_ms,dijkstra_heap_ms\n";
  std::cout << " n m Roy-Floyd n x Dijkstra, no heap"
            << " n x Dijkstra, heap\n";
  double prev[3] = {0, 0, 0};
  long long n = n0;
  for (int i = 0; i < count; i++, n *= 2) {
    long long m = std::max(n - 1, degree > 0 ? n * degree : n * (n - 1) / 4);
    Graph g;
    std::string err = "n too large";
    uint32_t seed = static_cast<uint32_t>(1 + i);
    bool ok = n <= MAX_N;
    if (!ok || !generate(g, static_cast<int>(n), m, 1, 100, seed, false, err)) {
      std::cout << " n = " << n << ": stopped (" << err << ")\n";
      break;
    }
    std::cout << std::setw(9) << n << std::setw(11) << m << std::flush;
    double t[3];
    for (int k = 0; k < 3; k++) {
      t[k] = timeMs(ALGS[k], g, nullptr);
      std::ostringstream cell;
      cell << timeStr(t[k]);
      if (prev[k] > 0)
        cell << " x" << std::fixed << std::setprecision(1) << t[k] / prev[k];
      std::cout << std::setw(k == 0   ? 18
                             : k == 1 ? 24
                                      : 20)
                << cell.str() << std::flush;
      prev[k] = t[k];
    }
    std::cout << '\n';
    csv << n << ',' << m << ',' << t[0] << ',' << t[1] << ',' << t[2] << '\n';
  }
  std::cout << " xk = T(2n) / T(n). The results are also in the file "
               "dublare_rf.csv.\n";
}

// ---------------------------------------------------------------- menu

enum Ask { OK, BAD, END }; // valid answer, invalid answer, end of input

static Ask askLine(const std::string &prompt, std::string &line) {
  std::cout << prompt;
  if (!std::getline(std::cin, line))
    return END;
  size_t a = line.find_first_not_of(" \t\r"),
         b = line.find_last_not_of(" \t\r");
  line = a == std::string::npos ? "" : line.substr(a, b - a + 1);
  return OK;
}

// exactly one integer in [lo, hi]; *dflt is used for an empty line
static Ask askNum(const std::string &prompt, long long lo, long long hi,
                  long long &x, const long long *dflt = nullptr) {
  std::string line;
  if (askLine(prompt, line) == END)
    return END;
  if (line.empty() && dflt) {
    x = *dflt;
    return OK;
  }
  std::vector<long long> v;
  if (!parseNums(line, v) || v.size() != 1 || v[0] < lo || v[0] > hi) {
    std::cout << "Invalid value: enter an integer between " << lo << " and "
              << hi << ".\n";
    return BAD;
  }
  x = v[0];
  return OK;
}

static Ask askYes(const std::string &prompt, bool &yes) {
  std::string line;
  if (askLine(prompt + " (y/n): ", line) == END)
    return END;
  if (line == "y" || line == "yes") {
    yes = true;
    return OK;
  }
  if (line == "n" || line == "no") {
    yes = false;
    return OK;
  }
  std::cout << "Answer y or n.\n";
  return BAD;
}

struct Session {
  Graph g;
  Apsp last; // the result of the last run of option 5
  bool ran = false;
};

static Ask optLoadSave(Session &S, bool load) {
  std::string file, err;
  if (askLine("File name: ", file) == END)
    return END;
  if (file.empty()) {
    std::cout << "Empty name.\n";
    return BAD;
  }
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
  S.ran = false;
  std::cout << "Loaded the graph.\n";
  printGraph(S.g);
  return OK;
}

static Ask optGenerate(Session &S) {
  long long n = 0, m = 0, a = 0, b = 0, seed = 0, one = 1;
  bool avoid = false;
  Ask r = askNum("Number of vertices n: ", 1, MAX_N, n);
  if (r == OK)
    r = askNum("Number of arcs m: ", 0, 5000000, m);
  if (r == OK)
    r = askNum("Minimum cost a: ", -1000000, 1000000, a);
  if (r == OK)
    r = askNum("Maximum cost b: ", -1000000, 1000000, b);
  if (r == OK)
    r = askNum("Seed (Enter = 1): ", 0, UINT32_MAX, seed, &one);
  if (r == OK && a < 0)
    r = askYes("No negative cycles?", avoid);
  if (r != OK)
    return r;
  Graph h;
  std::string err;
  if (!generate(h, static_cast<int>(n), m, static_cast<int>(a),
                static_cast<int>(b), static_cast<uint32_t>(seed), avoid, err)) {
    std::cout << "Invalid parameters: " << err << ".\n";
    return BAD;
  }
  S.g = h;
  S.ran = false;
  std::cout << "Generated a " << S.g.note << ".\n";
  return OK;
}

static void runRoyFloyd(Session &S, bool quiet) {
  double t = timeMs(royFloyd, S.g, &S.last);
  S.ran = true;
  std::cout << "Roy-Floyd: " << timeStr(t) << ", "
            << plural(S.last.compared, "comparison", "comparisons") << ", "
            << plural(S.last.updates, "value decreased", "values decreased")
            << ".\n";
  if (!quiet)
    printResult(S.g, S.last);
}

static void printPath(const Session &S, int s, int t) {
  const Graph &g = S.g;
  const Apsp &r = S.last;
  std::cout << " from " << g.label(s) << " to " << g.label(t) << ": ";
  if (r.d[s][t] == INF) {
    std::cout << "there is no path.\n";
    return;
  }
  if (r.d[s][t] == NEG_INF) {
    std::cout
        << "distance -inf: a path can go around a negative cycle any number "
        << "of times.\n";
    return;
  }
  std::vector<int> p = pathP(r, s, t), q = pathSplit(r, s, t);
  std::cout << "distance " << r.d[s][t] << ", path " << pathStr(g, p) << " ("
            << plural(static_cast<long long>(p.size()) - 1, "arc", "arcs")
            << ")\n";
  if (q.empty())
    std::cout
        << " the method without P failed: the graph has a cycle of cost 0\n";
  else if (q != p)
    std::cout << " the method without P finds another path, of the same cost: "
              << pathStr(g, q) << '\n';
}

static Ask optPath(Session &S) {
  const Graph &g = S.g;
  if (!S.ran)
    runRoyFloyd(S, true);
  long long x = 0, y = 0, none = -1;
  if (!g.queries.empty()) {
    Ask r = askNum("Start vertex (Enter = the pairs of the file): ", g.base,
                   g.base + g.n - 1, x, &none);
    if (r != OK)
      return r;
    if (x == none) {
      for (const auto &q : g.queries)
        printPath(S, q.first, q.second);
      return OK;
    }
  } else {
    Ask r = askNum("Start vertex: ", g.base, g.base + g.n - 1, x);
    if (r != OK)
      return r;
  }
  Ask r = askNum("End vertex: ", g.base, g.base + g.n - 1, y);
  if (r != OK)
    return r;
  printPath(S, static_cast<int>(x) - g.base, static_cast<int>(y) - g.base);
  return OK;
}

static void optWarshall(const Session &S) {
  const Graph &g = S.g;
  long long updates = 0;
  typedef std::chrono::steady_clock Clock;
  Clock::time_point t0 = Clock::now();
  Matrix<char> t = warshall(g, updates);
  double ms =
      std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
  long long ones = 0;
  for (char c : t.a)
    ones += c;
  std::cout << "Warshall: " << timeStr(ms) << ", "
            << plural(updates, "new pair", "new pairs") << "; " << ones - g.n
            << " of " << static_cast<long long>(g.n) * (g.n - 1)
            << " pairs (i, j), i != j, have j reachable from i.\n";
  if (g.n <= 40) {
    std::vector<std::string> cell;
    for (char c : t.a)
      cell.push_back(c ? "1" : "0");
    std::cout
        << " The transitive closure T (T[i][j] = 1: j is reachable from i):\n";
    printMatrix(g, cell, nullptr);
  }
}

static Ask optDoubling() {
  long long n0 = 0, count = 0, degree = 0, d0 = 100, c0 = 5, g0 = 5;
  Ask r = askNum("Initial number of vertices n0 (Enter = 100): ", 2, MAX_N, n0,
                 &d0);
  if (r == OK)
    r = askNum("Number of graphs (Enter = 5): ", 1, 10, count, &c0);
  if (r == OK)
    r = askNum("Mean degree m / n, 0 for a dense graph (Enter = 5): ", 0, 1000,
               degree, &g0);
  if (r != OK)
    return r;
  doubling(static_cast<int>(n0), static_cast<int>(count),
           static_cast<int>(degree));
  return OK;
}

int main(int argc, char **argv) {
  Session S;
  std::string err;
  std::cout << "All-pairs shortest paths, Roy-Floyd - DSA lab\n";
  if (argc > 1) {
    if (loadGraph(argv[1], S.g, err))
      std::cout << "Loaded " << argv[1] << ".\n";
    else
      std::cout << argv[1] << ": " << err << ".\n";
  }
  for (;;) {
    std::cout
        << "\n=== graph: " << plural(S.g.n, "vertex", "vertices") << ", "
        << plural(S.g.m(), "arc", "arcs") << " ===\n"
        << " 1. Load a graph from a file\n"
        << " 2. Generate a random graph\n"
        << " 3. Save the graph to a file\n"
        << " 4. Show the graph and the cost matrix\n"
        << " 5. Roy-Floyd: the distance matrix D and the matrix P, O(n^3)\n"
        << " 6. Roy-Floyd step by step: the matrix D after every vertex k\n"
        << " 7. The shortest path between two vertices\n"
        << " 8. Warshall: the transitive closure\n"
        << " 9. Compare Roy-Floyd with n runs of Dijkstra\n"
        << "10. The doubling experiment, on random graphs\n"
        << " 0. Exit\n";
    long long o = 0;
    Ask r = askNum("Option: ", 0, 10, o);
    if (r == END)
      break;
    if (r == BAD)
      continue;
    if (o == 0)
      break;
    if (o >= 3 && o <= 9 && S.g.n == 0) {
      std::cout << "There is no graph yet: load one (1) or generate one (2).\n";
      continue;
    }
    if (o == 1 || o == 3) {
      r = optLoadSave(S, o == 1);
    } else if (o == 2) {
      r = optGenerate(S);
    } else if (o == 4) {
      printGraph(S.g);
    } else if (o == 5) {
      runRoyFloyd(S, false);
    } else if (o == 6) {
      if (S.g.n > 12)
        std::cout << "The steps are shown for at most 12 vertices.\n";
      else
        traceRoyFloyd(S.g);
    } else if (o == 7) {
      r = optPath(S);
    } else if (o == 8) {
      optWarshall(S);
    } else if (o == 9) {
      if (S.g.hasNegative())
        std::cout
            << "Note: the graph has negative-cost arcs, so Dijkstra may give "
            << "wrong distances.\n";
      compareAll(S.g);
    } else {
      r = optDoubling();
    }
    if (r == END)
      break;
  }
  std::cout << "Goodbye!\n";
  return 0;
}