// Hungarian Algorithm, O(n^2 m) 

#include <bits/stdc++.h>

using namespace std;

typedef long long ll;

const ll INF = 1e16 + 5;

pair <ll, vector <int>> hungarian (const vector <vector <ll>> &a) {
  if (a.empty()) return {0, {}};
  int n = a.size() + 1, m = a[0].size() + 1;
  vector <ll> u(n), v(m);
  vector <int> p(m), ans(n - 1);
  for (int i = 1; i < n; ++i) {
    p[0] = i; int x = 0, y, z;
    vector <int> pre(m, -1);
    vector <bool> done(m + 1);
    vector <ll> dist(m, INF);
    do {
      done[x] = 1;
      y = p[x], delta = INF;
      for (int j = 1; j < m; ++j) if (!done[j]) {
        ll cur = a[y - 1][j - 1] - u[y] - v[j];
        if (cur < dist[j]) dist[j] = cur, pre[j] = x;
        if (dist[j] < delta) delta = dist[j], z = j;
      }
      for (int j = 0; j < m; ++j) {
        if (done[j]) u[p[j]] += delta, v[j] -= delta;
        else dist[j] -= delta;
      } x = z;
    } while (p[x]);
    while (x) z = pre[x], p[x] = p[z], x = z;
  }
  for (int j = 1; j < m; ++j) if (p[j]) ans[p[j] - 1] = j - 1;
  return {-v[0], ans};
}

int main() {

  return 0;
}

