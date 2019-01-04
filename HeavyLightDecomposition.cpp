// Solves UVa 12424

#include <bits/stdc++.h>

using namespace std;

const int C = 12;
const int N = 100010;

struct FenwickTree {
  int f[N];

  void clear() {
    memset(f, 0, sizeof f);
  }

  void update (int p, int v) {
    while (p < N) f[p] += v, p += p & -p;
  }

  int query (int p) {
    int res = 0;
    while (p > 0) res += f[p], p -= p & -p;
    return res;
  }

  int query (int l, int r) {
    return l > r ? 0 : query(r) - query(l - 1);
  }
};

FenwickTree f[C];
vector <int> g[N];
int t, n, m, ptr, c[N], par[N];
int sz[N], h[N], in[N], nxt[N];

void dfs (int u = 1, int far = 0) {
  sz[u] = 1, h[u] = far;
  for (int v : g[u]) g[v].erase(find(g[v].begin(), g[v].end(), u));
  for (int &v : g[u]) {
    par[v] = u; dfs(v, far + 1); sz[u] += sz[v];
    if (sz[v] > sz[g[u][0]]) swap(v, g[u][0]);
  }
}

void hld (int u = 1) {
  in[u] = ++ptr;
  for (int v : g[u]) {
    nxt[v] = (v == g[u][0] ? nxt[u] : v); hld(v);
  }
}

int query (int col, int u, int v) {
  int res = 0;
  while (nxt[u] != nxt[v]) {
    if (h[nxt[u]] > h[nxt[v]]) swap(u, v);
    res += f[col].query(in[nxt[v]], in[v]);
    v = par[nxt[v]];
  }
  if (h[u] > h[v]) swap(u, v);
  res += f[col].query(in[u], in[v]);
  return res;
}

int main() {
  cin >> t;
  while (t--) {
    cin >> n >> m;

    for (int i = 1; i <= n; ++i) g[i].clear();
    for (int i = 1; i < C; ++i) f[i].clear();
    
    for (int i = 1; i <= n; ++i) scanf("%d", c + i);
    for (int i = 1; i < n; ++i) {
      int u, v; scanf("%d %d", &u, &v);
      g[u].push_back(v), g[v].push_back(u);
    }
    
    ptr = 0; dfs(); hld();
    for (int i = 1; i <= n; ++i) f[c[i]].update(in[i], +1);
    
    while (m--) {
      int cmd, u, v;
      scanf("%d %d %d", &cmd, &u, &v);
      if (cmd) {
        int ans = 0;
        for (int i = 1; i < C; ++i) ans = max(ans, query(i, u, v));
        printf("%d\n", ans);
      } else {
        f[c[u]].update(in[u], -1);
        c[u] = v;
        f[c[u]].update(in[u], +1);
      }
    }
  }
  return 0;
}

