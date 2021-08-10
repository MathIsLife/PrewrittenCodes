// O(n^2 2^n) subset convolution
// ans[i] = sum f[j] g[i ^ j] over subsets j of i

#include <bits/stdc++.h>

using namespace std;

const int LG = 18;
const int N = 1 << LG;

int n, f[N], g[N], f_hat[LG][N], g_hat[LG][N], h[LG][N], ans[N];

int main() {
  for (int mask = 0; mask < 1 << n; ++mask) {
    f_hat[__builtin_popcount(mask)][mask] = f[mask];
    g_hat[__builtin_popcount(mask)][mask] = g[mask];
  }
  for (int i = 0; i <= n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int mask = 0; mask < 1 << n; ++mask) {
        if (mask & 1 << j) {
          f_hat[i][mask] += f_hat[i][mask ^ 1 << j];
          g_hat[i][mask] += g_hat[i][mask ^ 1 << j];
        }
      }
    }
  }
  for (int mask = 0; mask < 1 << n; ++mask) {
    for (int i = 0; i <= n; ++i) {
      for (int j = 0; j <= i; ++j) {
        h[i][mask] += f_hat[j][mask] * g_hat[i - j][mask];
      } 
    }
  }
  for (int i = 0; i <= n; ++i) {
    for (int j = 0; j < n; ++j) {
      for (int mask = 0; mask < 1 << n; ++mask) {
        if (mask & 1 << j) {
          h[i][mask] -= h[i][mask ^ 1 << j];
        }
      }
    }
  }
  for (int mask = 0; mask < 1 << n; ++mask) ans[mask] = h[__builtin_popcount(mask)][mask];
  return 0;
}

