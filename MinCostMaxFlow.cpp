#include <bits/stdc++.h>

using namespace std;

namespace MCMF {
  using T = long long;

  const int MAX = 100010;
  const T INF = 1e12 + 69;

  bool vis[MAX];
  T dist[MAX], caps[MAX];
  int n, par[MAX], pos[MAX];

  struct edge {
    int to, rev_pos;
    T cap, cost, flow;
  };

  vector <edge> g[MAX];

  inline void init (int nodes) {
    n = nodes;
    for (int i = 0; i < n; ++i) g[i].clear();
  }

  inline void AddEdge (int u, int v, T cap, T cost) {
    edge a = {v, g[v].size(), cap, cost, 0};
    edge b = {u, g[u].size(), 0, -cost, 0};
    g[u].emplace_back(a);
    g[v].emplace_back(b);
  }

  bool SPFA (int src, int snk) {
    for (int i = 0; i < n; ++i) {
      caps[i] = dist[i] = INF, vis[i] = 0;
    }
    queue <int> q;
    dist[src] = 0, vis[src] = 1, q.emplace(src);
    while (!q.empty()) {
      int u = q.front(); q.pop();
      vis[u] = 0;
      for(int i = 0; i < g[u].size(); ++i) {
        edge &e = g[u][i];
        int v = e.to;
        if (e.cap > e.flow and dist[v] > dist[u] + e.cost){
          dist[v] = dist[u] + e.cost, par[v] = u, pos[v] = i;
          caps[v] = min(caps[u], e.cap - e.flow);
          if (!vis[v]) vis[v] = 1, q.emplace(v);
        }
      }
    }
    return dist[snk] != INF;
  }

  pair <T, T> MinCostFlow (int src, int snk) {
    int u, v;
    T flow = 0, cost = 0, f;
    while (SPFA(src, snk)) {
      u = snk, f = caps[u], flow += f;
      while (u ^ src){
        v = par[u];
        g[v][pos[u]].flow += f;
        g[u][g[v][pos[u]].rev_pos].flow -= f;
        u = v;
      }
      cost += dist[snk] * f;
    }
    return make_pair(flow, cost);
  }
}

int main() {
  int n, m;
  cin >> n >> m;
  MCMF::init(n);
  while (m--) {
    int u, v, cap, cost;
    cin >> u >> v >> cap >> cost;
    MCMF::AddEdge(u, v, cap, cost);
  }
  auto [flow, cost] = MCMF::MinCostFlow(0, n - 1);
  cout << flow << " " << cost << '\n';
  return 0;
}

