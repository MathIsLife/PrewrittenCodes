#include <bits/stdc++.h>

using namespace std;

const int N = 505;

struct Blossom {
  queue <int> q;
  vector <int> g[N];
  int t, n, vis[N], par[N], orig[N], match[N], aux[N];

  Blossom (int _n = 0) {
    n = _n, t = 0;
    for (int i = 0; i <= n; ++i) {
      g[i].clear(), match[i] = aux[i] = par[i] = 0;
    }
  }

  inline void addEdge (int u, int v) {
    g[u].push_back(v), g[v].push_back(u);
  }

  void augment (int u, int v) {
    int pv = v, nv;
    do {
      pv = par[v], nv = match[pv];
      match[v] = pv, match[pv] = v, v = nv;
    } while (u ^ pv);
  }

  int lca (int u, int v) {
    ++t;
    while (true) {
      if (u) {
        if (aux[u] == t) return u; aux[u] = t;
        u = orig[par[match[u]]];
      }
      swap(u, v);
    }
  }

  void blossom (int u, int v, int x) {
    while (orig[u] ^ x) {
      par[u] = v, v = match[u];
      if (vis[v] == 1) q.emplace(v), vis[v] = 0;
      orig[u] = orig[v] = x, u = par[v];
    }
  }

  bool bfs (int src) {
    fill(vis + 1, vis + n + 1, -1); 
    iota(orig + 1, orig + n + 1, 1);
    while (!q.empty()) q.pop();
    q.emplace(src), vis[src] = 0;
    while (!q.empty()) {
      int u = q.front(); q.pop();
      for (int v : g[u]) {
        if (vis[v] == -1) {
          par[v] = u, vis[v] = 1;
          if (!match[v]) return augment(src, v), 1;
          q.emplace(match[v]), vis[match[v]] = 0;
        } else if (vis[v] == 0 and orig[u] ^ orig[v]) {
          int x = lca(orig[u], orig[v]);
          blossom(v, u, x), blossom(u, v, x);
        }
      }
    } return 0;
  }

  int maxMatch() {
    int ans = 0;
    vector <int> vec(n - 1); 
    iota(vec.begin(), vec.end(), 1);
    shuffle(vec.begin(), vec.end(), mt19937(69));
    for (int u : vec) if (!match[u]) {
      for (int v : g[u]) if (!match[v]) {
        match[u] = v, match[v] = u;
        ++ans; break;
      }
    }
    for (int i = 1; i <= n; ++i) if (!match[i] and bfs(i)) ++ans;
    return ans;
  }
};

int n, m;

int main() {
  cin >> n >> m;
  Blossom yo(n);
  while (m--) {
    int u, v;
    scanf("%d %d", &u, &v);
    ++u, ++v;
    yo.addEdge(u, v);
  }
  cout << yo.maxMatch() << '\n';
  for (int i = 1; i <= n; ++i) {
    if (yo.match[i] > i) printf("%d %d\n", i - 1, yo.match[i] - 1);
  }
  return 0;
}

