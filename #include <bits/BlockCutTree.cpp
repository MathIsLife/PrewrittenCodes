#include <bits/stdc++.h>

using namespace std;

const int N = 200010;

bitset <N> art, good;
vector <int> g[N], tree[N], st, comp[N];
int n, m, ptr, cur, in[N], low[N], id[N];

void dfs (int u, int from = -1) {
  in[u] = low[u] = ++ptr;
  st.emplace_back(u);
  for (int v : g[u]) if (v ^ from) {
    if (!in[v]) {
      dfs(v, u);
      low[u] = min(low[u], low[v]);
      if (low[v] >= in[u]) {
        art[u] = in[u] > 1 or in[v] > 2;
        comp[++cur].emplace_back(u);
        while (comp[cur].back() ^ v) {
          comp[cur].emplace_back(st.back());
          st.pop_back();
        }
      } 
    } else {
      low[u] = min(low[u], in[v]);
    }
  }
}

void buildTree() {
  ptr = 0;
  for (int i = 1; i <= n; ++i) {
    if (art[i]) id[i] = ++ptr;
  }
  for (int i = 1; i <= cur; ++i) {
    int x = ++ptr;
    for (int u : comp[i]) {
      if (art[u]) {
        tree[x].emplace_back(id[u]);
        tree[id[u]].emplace_back(x);
      } else {
        id[u] = x;
      }
    }
  }
}

int main() {
  cin >> n >> m;
  while (m--) {
    int u, v;
    scanf("%d %d", &u, &v);
    g[u].emplace_back(v);
    g[v].emplace_back(u);
  }
  for (int i = 1; i <= n; ++i) {
    if (!in[i]) dfs(i);
  }
  buildTree();
  for (int i = 1; i <= ptr; ++i) {
    cout << i << " --> ";
    for (int j : tree[i]) cout << j << " "; cout << '\n';
  }
  return 0;
}

