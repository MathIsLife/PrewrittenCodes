// Solves SPOJ FASTFLOW

#include <bits/stdc++.h>

using namespace std;

typedef long long ll;

struct edge {
  int dest, back; ll f, c;

  edge (int d, int b, ll f, ll c) : dest(d), back(b), f(f), c(c) {}
};

struct PushRelabel {
  vector <vector <edge>> g;
  vector <ll> ec;
  vector <edge*> cur;
  vector <vector <int>> hs; 
  vector <int> H;
  
  PushRelabel (int n) : g(n), ec(n), cur(n), hs(n + n), H(n) {}

  void AddEdge (int s, int t, ll cap, ll rcap = 0) {
    if (s == t) return;
    g[s].emplace_back(t, g[t].size(), 0, cap);
    g[t].emplace_back(s, g[s].size() - 1, 0, rcap);
  }

  void AddFlow (edge &e, ll f) {
    edge &back = g[e.dest][e.back];
    if (!ec[e.dest] and f) hs[H[e.dest]].emplace_back(e.dest);
    e.f += f; e.c -= f; ec[e.dest] += f;
    back.f -= f; back.c += f; ec[back.dest] -= f;
  }

  ll MaxFlow (int s, int t) {
    int v = g.size(); H[s] = v; ec[t] = 1;
    vector <int> co(v + v); co[0] = v - 1;
    for (int i = 0; i < v; ++i) cur[i] = g[i].data();
    for (auto &e : g[s]) AddFlow(e, e.c);
    for (int hi = 0; ;) {
      while (hs[hi].empty()) if (!hi--) return -ec[s];
      int u = hs[hi].back(); hs[hi].pop_back();
      while (ec[u] > 0) 
        if (cur[u] == g[u].data() + g[u].size()) {
          H[u] = 1e9 + 5;
          for (auto &e : g[u]) if (e.c and H[u] > H[e.dest] + 1)
            H[u] = H[e.dest] + 1, cur[u] = &e;
          if (++co[H[u]], !--co[hi] and hi < v)
            for (int i = 0; i < v; ++i) if (hi < H[i] and H[i] < v)
              --co[H[i]], H[i] = v + 1;
          hi = H[u];
        } else if (cur[u] -> c and H[u] == H[cur[u] -> dest] + 1)
          AddFlow(*cur[u], min(ec[u], cur[u] -> c));
        else ++cur[u];
    }
  }
};

int main() {
  int n, m;
  cin >> n >> m;
  PushRelabel flow(n);
  while (m--) {
    int u, v, c;
    scanf("%d %d %d", &u, &v, &c);
    --u, --v;
    flow.AddEdge(u, v, c, c);
  }
  cout << flow.MaxFlow(0, n - 1) << '\n';
  return 0;
}

