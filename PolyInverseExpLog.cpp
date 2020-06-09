#include <bits/stdc++.h>

using namespace std;

typedef long long ll;

const int G = 3;
const int MOD = 998244353;
const int N = (1 << 20) + 5;

int rev[N], w[N], inv_n, inv[N];

int bigMod (int a, int e) {
  if (e == -1) e = MOD - 2;
  int ret = 1;
  while (e) {
    if (e & 1) ret = (ll) ret * a % MOD;
    a = (ll) a * a % MOD; e >>= 1;
  }
  return ret;
}

void prepare (int n) { 
  int sz = abs(31 - __builtin_clz(n));
  int r = bigMod(G, (MOD - 1) / n); 
  inv_n = inv[n], w[0] = w[n] = 1; 
  for (int i = 1; i < n; ++i) w[i] = (ll) w[i - 1] * r % MOD;
  for (int i = 1; i < n; ++i) rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (sz - 1));
}

static int f[N], g[N], h[N];

void ntt (int *a, int n, int dir = 0) { 
  for (int i = 1; i < n - 1; ++i) { 
    if (i < rev[i]) swap(a[i], a[rev[i]]);
  }
  for (int m = 2; m <= n; m <<= 1) {
    for (int i = 0; i < n; i += m) {
      for (int j = 0; j < (m >> 1); ++j) {
        int &u = a[i + j], &v = a[i + j + (m >> 1)]; 
        int t = (ll) v * w[dir ? n - n / m * j : n / m * j] % MOD;
        v = u - t < 0 ? u - t + MOD : u - t;
        u = u + t >= MOD ? u + t - MOD : u + t;
      }
    }
  } if (dir) for (int i = 0; i < n; ++i) a[i] = (ll) a[i] * inv_n % MOD;
}

void multiply (int *a, int *b, int n) {
  prepare(n << 1), ntt(a, n << 1), ntt(b, n << 1);
  for (int i = 0; i < n << 1; ++i) a[i] = (ll) a[i] * b[i] % MOD;
  ntt(a, n << 1, 1); for (int i = n; i < n << 1; ++i) a[i] = 0;
}

void inverse (int *a, int n, int *b) {
  b[0] = inv[a[0]], b[1] = 0;
  for (int m = 2; m <= n; m <<= 1) {
    for (int i = 0; i < m; ++i) f[i] = a[i], f[i + m] = b[i + m] = 0;
    prepare(m << 1), ntt(f, m << 1), ntt(b, m << 1);
    for (int i = 0; i < m << 1; ++i) b[i] = (ll) b[i] * (MOD + 2 - (ll) b[i] * f[i] % MOD) % MOD;
    ntt(b, m << 1, 1); for (int i = m; i < m << 1; ++i) b[i] = 0;
  }
}

void log (int *a, int n) {
  inverse(a, n, g);
  for (int i = 0; i + 1 < n; ++i) a[i] = (i + 1LL) * a[i + 1] % MOD;
  multiply(a, g, n); 
  for (int i = n - 1; i; --i) a[i] = (ll) a[i - 1] * inv[i] % MOD;
  a[0] = 0;
}

void exp (int *a, int n, int *b) {
  b[0] = 1, b[1] = 0;
  for (int m = 2; m <= n; m <<= 1) {
    for (int i = 0; i < m; ++i) h[i] = b[i];
    log(h, m);
    for (int i = 0; i < m; ++i) h[i] = (a[i] - h[i] + MOD) % MOD;
    ++h[0], h[0] %= MOD;
    for (int i = m; i < m << 1; ++i) b[i] = h[i] = 0;
    multiply(b, h, m); for (int i = m; i < m << 1; ++i) b[i] = 0;
  }
}

int n, a[N], b[N];

int main() {
  cin >> n;
  for (int i = 0; i < n; ++i) scanf("%d", a + i);
  int m = 1; while (m < n) m <<= 1;
  exp(a, m, b);
  for (int i = 0; i < n; ++i) printf("%d ", b[i]); puts("");
  return 0;
}

