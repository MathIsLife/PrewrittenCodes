// Shortest cycle in undirected graph: O(VE)
// Solves LA 6151 - Beehives

#include <bits/stdc++.h>
 
using namespace std;
 
const int N = 505;
const int M = 5000010;
const int INF = 1e9 + 5;
 
bitset <N> done;
vector <int> g[N];
int t, cs, n, m, d[N], ans, par[N], q[M];
 
void relax (int src) {
  done.reset();
  par[src] = -1, d[src] = 0;
  int st = 0, en = 0; q[0] = src;
  while (st <= en) {
    int u = q[st++]; done[u] = 1;
    if (d[u] > (ans + 1) / 2) return;
    for (int v : g[u]) if (v - par[u]) {
      if (!done[v]) {
        d[v] = d[u] + 1, par[v] = u;
        if (d[v] < ans) q[++en] = v;
      } else {
        ans = min(ans, d[u] + d[v] + 1);
      }
    }
  }
}
 
int main() {
  cin >> t;
  while (t--) {
    scanf("%d %d", &n, &m);
    while (m--) {
      int u, v;
      scanf("%d %d", &u, &v);
      g[u].emplace_back(v);
      g[v].emplace_back(u);
    }
    ans = INF;
    for (int i = 0; i < n; ++i) relax(i);
    printf("Case %d: ", ++cs);
    if (ans == INF) puts("impossible");
    else printf("%d\n", ans);
    for (int i = 0; i < n; ++i) g[i].clear();
  }
  return 0;
}

