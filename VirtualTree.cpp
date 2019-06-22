#include <bits/stdc++.h>

using namespace std;

const int LG = 19;
const int N = 100010;

vector <int> g[N], virt[N], cost[N];
int n, m, ptr, h[N], in[N], p[N][LG], stk[N];

void go (int u = 1, int from = -1, int far = 0) {
  h[u] = far, p[u][0] = from, in[u] = ++ptr;
  for (int v : g[u]) if (v - from) go(v, u, far + 1);
}

int get_lca (int u, int v) {
  if (h[u] < h[v]) swap(u, v);
  for (int i = LG - 1; i >= 0; --i) {
    if (h[u] - (1 << i) >= h[v]) u = p[u][i];
  }
  if (u == v) return u;
  for (int i = LG - 1; i >= 0; --i) {
    if (~p[u][i] and p[u][i] - p[v][i]) {
      u = p[u][i], v = p[v][i];
    }
  }
  return p[u][0];
}

void add_edge (int u, int v) {
  if (u == v) return;
  virt[u].emplace_back(v);
  virt[v].emplace_back(u);
  int w = abs(h[u] - h[v]);
  cost[u].emplace_back(w);
  cost[v].emplace_back(w);
}

void buildTree (vector <int> &nodes) {
  if (nodes.size() <= 1) return;
  sort(nodes.begin(), nodes.end(), [] (int x, int y) {return in[x] < in[y];});
  int root = get_lca(nodes[0], nodes.back()), sz = nodes.size();
  ptr = 0, stk[ptr++] = root;
  for (int i = 0; i < sz; ++i) {
    int u = nodes[i], lca = get_lca(u, stk[ptr - 1]);
    if (lca == stk[ptr - 1]) {
      stk[ptr++] = u;
    } else {
      while (ptr > 1 and h[stk[ptr - 2]] >= h[lca]) {
        add_edge(stk[ptr - 2], stk[ptr - 1]), --ptr;
      }
      if (stk[ptr - 1] != lca) {
        add_edge(lca, stk[--ptr]);
        stk[ptr++] = lca, nodes.emplace_back(lca);
      }
      stk[ptr++] = u;
    }
  }
  if (find(nodes.begin(), nodes.end(), root) == nodes.end()) nodes.emplace_back(root);
  for (int j = 0; j + 1 < ptr; ++j) add_edge(stk[j], stk[j + 1]);
}

int main() {
  cin >> n;
  for (int i = 1, u, v; i < n; ++i) {
    scanf("%d %d", &u, &v);
    g[u].emplace_back(v);
    g[v].emplace_back(u);
  }
  memset(p, -1, sizeof p); go();
  for (int j = 1; j < LG; ++j) {
    for (int i = 1; i <= n; ++i) {
      if (~p[i][j - 1]) p[i][j] = p[p[i][j - 1]][j - 1];
    }
  }
  cin >> m;
  vector <int> nodes(m);
  for (int i = 0; i < m; ++i) {
    scanf("%d", &nodes[i]);
  }
  buildTree(nodes);
  for (int u : nodes) {
    cout << u << " --> ";
    for (int v : virt[u]) cout << v << " ";
    cout << '\n';
  }
  return 0;
}

