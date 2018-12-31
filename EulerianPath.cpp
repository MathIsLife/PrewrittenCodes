#include <bits/stdc++.h>

using namespace std;

// Eulerian path / circuit

// Undirected graph: circuit (or edge disjoint cycles) exists iff all nodes are of even degree
// Undirected graph: path exists iff number of odd degree nodes is zero or two

// Directed graph: circuit (or edge disjoint directed cycles) exists iff each node 
//   satisfies in_degree = out_degree and the graph is strongly connected
// Directed graph: path exists iff at most one vertex has in_degree - out_degree = 1
//   and at most one vertex has out_degree - in_degree = 1 and all other vertices have
//   in_degree = out_degree, and graph is weakly connected

const int N = 200010;

bitset <N> bad;
vector <int> g[N];
vector <int> circ;
int n, m, deg[N], U[N], V[N];

void hierholzer (int src) {
  if (!deg[src]) return;
  vector <int> path;
  path.push_back(src);
  int at = src;
  while (!path.empty()) {
    if (deg[at]) {
      path.push_back(at);
      while (bad[g[at].back()]) g[at].pop_back();
      int e = g[at].back(), nxt = U[e] ^ at ^ V[e];
      bad[e] = 1, --deg[at], --deg[nxt], at = nxt;
    } else {
      circ.push_back(at);
      at = path.back(), path.pop_back();
    }
  }
  reverse(circ.begin(), circ.end());
}

int main() {
  cin >> n >> m;
  for (int i = 1; i <= m; ++i) {
    scanf("%d %d", U + i, V + i);
    g[U[i]].push_back(i); 
    g[V[i]].push_back(i);
    ++deg[U[i]], ++deg[V[i]];
  }
  hierholzer(1);
  for (int x : circ) printf("%d ", x); puts("");
  return 0;
}

