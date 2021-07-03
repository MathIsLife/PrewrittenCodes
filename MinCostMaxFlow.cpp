#include <bits/stdc++.h>

using namespace std;

typedef long long ll;

const int N = 1234;
const ll INF = 1e9 + 69;

struct edge {
  int v, rev_pos;
  ll cap, cost, flow = 0;
};

struct MCMF {
  int n;
  vector <vector <edge>> E;

  bitset <N> vis;
  int par[N], pos[N];
  ll caps[N], dist[N];

  MCMF (int n) : n(n), E(n) {}

  inline void AddEdge (int u, int v, ll cap, ll cost = 0) {
    E[u].push_back({v, E[v].size(), cap, cost});
    E[v].push_back({u, E[u].size() - 1, 0, -cost});
  }

  bool SPFA (int S, int T) {
    vis.reset();
    for(int i = 0; i < n; ++i) caps[i] = dist[i] = INF;
    queue <int> q;
    q.emplace(S), dist[S] = 0, vis[S] = 1;
    while (!q.empty()) {
      int i = 0, u = q.front();
      q.pop(), vis[u] = 0;
      for (auto &e : E[u]) {
        int v = e.v;
        ll f = e.cap - e.flow, w = dist[u] + e.cost;
        if (f > 0 and dist[v] > w) {
           dist[v] = w, par[v] = u, pos[v] = i;
           caps[v] = min(caps[u], f);
           if (!vis[v]) q.push(v), vis[v] = 1;
        } ++i;
      }
    }
    return dist[T] != INF;
  }

  pair <ll, ll> GetMaxFlow (int S, int T) {
    ll flow = 0, cost = 0;
    while (SPFA(S, T)) {
      int v = T;
      ll f = caps[v];
      flow += f;
      while (v ^ S) {
        int u = par[v];
        edge &e = E[u][pos[v]];
        e.flow += f;
        E[v][e.rev_pos].flow -= f;
        v = u;
      }
      cost += dist[T] * f;
    }
    return make_pair(flow, cost);
  }
};

int main() {
  int n, m;
  cin >> n >> m;
  MCMF flow_graph(n);
  while (m--) {
    int u, v, cap, cost;
    cin >> u >> v >> cap >> cost;
    flow_graph.AddEdge(u, v, cap, cost);
  }
  auto [flow, cost] = flow_graph.GetMaxFlow(0, n - 1);
  cout << flow << " " << cost << '\n';
  return 0;
}

