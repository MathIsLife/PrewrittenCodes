// Diameter with online positive edge weight updates, O(n lg n)

#include <bits/stdc++.h>
 
using namespace std;
 
typedef long long ll;

// f[x] = distance from root to x
// x <= y <= z
// max{f[x] - 2f[y] + f[z]}
// a --> f[x]
// b --> -2f[y]
// c --> f[x] - 2f[y]
// d --> -2f[y] + f[z]
// e --> f[x] - 2f[y] + f[z]
struct node {
  ll a = 0, b = 0, c = 0, d = 0, e = 0, lazy = 0;
 
  node operator + (const node &oth) const {
    node ret;
    ret.a = max(a, oth.a);
    ret.b = max(b, oth.b); 
    ret.c = max(max(c, oth.c), a + oth.b);
    ret.d = max(max(d, oth.d), b + oth.a);
    ret.e = max(max(e, oth.e), max(c + oth.a, a + oth.d));
    return ret;
  }
};
 
const int N = 200010;
 
ll mod, W[N];
node t[4 * N];
vector <int> g[N];
int n, q, U[N], V[N], ptr, l[N], r[N], pos[N];
 
void dfs (int u = 1, int from = 0) {
  l[u] = ++ptr;
  for (int e : g[u]) {
    int v = U[e] ^ u ^ V[e];
    if (v == from) continue;
    pos[e] = v; dfs(v, u); ++ptr;
  }
  r[u] = ptr;
}
 
void push (int u, int b, int e) {
  ll v = t[u].lazy;
  t[u].a += v, t[u].b -= v + v, t[u].c -= v, t[u].d -= v;
  if (b ^ e) t[u << 1].lazy += v, t[u << 1 | 1].lazy += v;
  t[u].lazy = 0;
}
 
void update (int l, int r, ll v, int u = 1, int b = 1, int e = n + n) {
  if (t[u].lazy) push(u, b, e);
  if (b > r or e < l) return;
  if (b >= l and e <= r) {
    t[u].lazy += v;
    return push(u, b, e);
  }
  int mid = b + e >> 1;
  update(l, r, v, u << 1, b, mid);
  update(l, r, v, u << 1 | 1, mid + 1, e);
  t[u] = t[u << 1] + t[u << 1 | 1];
}
 
int main() {
  cin >> n >> q >> mod;
  for (int i = 1; i < n; ++i) {
    scanf("%d %d %lld", U + i, V + i, W + i);
    g[U[i]].emplace_back(i), g[V[i]].emplace_back(i);
  }
  dfs();
  for (int i = 1; i < n; ++i) update(l[pos[i]], r[pos[i]], W[i]);
  ll last = 0;
  while (q--) {
    int d; ll e;
    scanf("%d %lld", &d, &e);
    d = 1 + (d + last) % (n - 1);
    e = (e + last) % mod;
    update(l[pos[d]], r[pos[d]], e - W[d]);
    last = t[1].e, W[d] = e;
    printf("%lld\n", last);
  }
  return 0;
}

