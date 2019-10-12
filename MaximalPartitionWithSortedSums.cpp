// Maximum number of partitions of an array so that the sums are sorted
// Array values are positive, O(n)

#include <bits/stdc++.h>

using namespace std;

const int N = 123456;

int n, at, ptr = 1, h[N], f[N], to[N], q[N];

int main() {
  cin >> n;
  for (int i = 1; i <= n; ++i) scanf("%d", h + i);
  for (int i = 1; i <= n; ++i) {
    h[i] += h[i - 1], f[i] = h[i] - h[q[at]], to[i] = q[at];
    while (at < ptr and f[q[at + 1]] + h[q[at + 1]] <= h[i]) {
      f[i] = h[i] - h[q[++at]], to[i] = q[at];
    }
    while (ptr and f[q[ptr]] + h[q[ptr]] > f[i] + h[i]) --ptr;
    q[++ptr] = i, at = min(at, ptr);
  }
  int ans = n;
  for (int i = n; i; i = to[i]) --ans;
  cout << ans << endl;
  return 0;
}

