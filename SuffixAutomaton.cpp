// Solves pattern matching problem (UVa 10679)

#include <bits/stdc++.h>

using namespace std;

const int N = 100010;

char s[N], p[N];
map <char, int> to[N << 1];
int len[N << 1], link[N << 1], sz, last;

inline void init() {
  len[0] = 0, link[0] = -1, sz = 1, last = 0, to[0].clear();
}

void feed (char c) {
  int cur = sz++, p = last;
  len[cur] = len[last] + 1, link[cur] = 0, to[cur].clear();
  while (~p and !to[p].count(c)) to[p][c] = cur, p = link[p];
  if (~p) {
    int q = to[p][c];
    if (len[q] - len[p] - 1) {
      int r = sz++;
      len[r] = len[p] + 1, to[r] = to[q], link[r] = link[q];
      while (~p and to[p][c] == q) to[p][c] = r, p = link[p];
      link[q] = link[cur] = r; 
    } else link[cur] = q;
  } last = cur;
}

bool run() {
  int m = strlen(p);
  for (int i = 0, u = 0; i < m; ++i) {
    if (!to[u].count(p[i])) return 0;
    u = to[u][p[i]];
  } return 1;
}

int main() {
  int t; cin >> t;
  while (t--) {
    scanf("%s", s);
    int n = strlen(s);
    init();
    for (int i = 0; i < n; ++i) feed(s[i]);
    int q; scanf("%d", &q);
    while (q--) {
      scanf("%s", p);
      puts(run() ? "y" : "n");
    }
  }
  return 0;
}

