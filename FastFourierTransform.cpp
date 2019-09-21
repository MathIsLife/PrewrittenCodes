#include <bits/stdc++.h>

using namespace std;

typedef long long ll;
typedef long double ld;

namespace FFT {
  const int T = 1 << 17;

  struct cplx {
    ld x, y;
    
    cplx (ld x = 0, ld y = 0) : x(x), y(y) {}

    cplx operator + (const cplx &rhs) const {return cplx(x + rhs.x, y + rhs.y);}
    cplx operator - (const cplx &rhs) const {return cplx(x - rhs.x, y - rhs.y);}
    cplx operator * (const cplx &rhs) const {return cplx(x * rhs.x - y * rhs.y, x * rhs.y + y * rhs.x);}
    cplx operator !() const {return cplx(x, -y);}
  };

  int rev[T];
  cplx rts[T + 1];
  cplx f_a[T], f_b[T];
  cplx f_c[T], f_d[T];

  void init() {
    int k = 0; while (1 << k < T) ++k;
    rev[0] = 0;
    for (int i = 1; i < T; ++i) {
      rev[i] = rev[i >> 1] >> 1 | ((i & 1) << k - 1);
    }
    ld PI = acosl(-1.0);
    rts[0] = rts[T] = cplx(1, 0);
    for (int i = 1; i + i <= T; ++i) {
      rts[i] = cplx(cos(i * 2 * PI / T), sin(i * 2 * PI / T));
    }
    for (int i = T / 2 + 1; i < T; ++i) {
      rts[i] = !rts[T - i];
    }
  }

  void dft (cplx a[], int n, int sign) {
    static int is_init;
    if (!is_init) is_init = 1, init();
    int d = 0; while ((1 << d) * n != T) ++d;
    for (int i = 0; i < n; ++i) {
      if (i < (rev[i] >> d)) {
        swap(a[i], a[rev[i] >> d]);
      }
    }
    for (int len = 2; len <= n; len <<= 1) {
      int delta = T / len * sign;
      for (int i = 0; i < n; i += len) {
        cplx *x = a + i, *y = a + i + (len >> 1), *w = sign > 0 ? rts : rts + T;
        for (int k = 0; k + k < len; ++k) {
          cplx z = *y * *w;
          *y = *x - z, *x = *x + z;
          ++x, ++y, w += delta;
        }
      }
    }
    if (sign < 0) {
      for (int i = 0; i < n; ++i) {
        a[i].x /= n, a[i].y /= n;
      }
    }
  }

  void multiply (int a[], int b[], int n_a, int n_b, ll c[], int dup = 0) {
    int n = n_a + n_b - 1; while (n != (n & -n)) n += n & -n;
    for (int i = 0; i < n; ++i) f_a[i] = f_b[i] = cplx();
    for (int i = 0; i < n_a; ++i) f_a[i] = cplx(a[i]);
    for (int i = 0; i < n_b; ++i) f_b[i] = cplx(b[i]);
    dft(f_a, n, 1);
    if (dup) for (int i = 0; i < n; ++i) f_b[i] = f_a[i];
    else dft(f_b, n, 1);
    for (int i = 0; i < n; ++i) f_a[i] = f_a[i] * f_b[i];
    dft(f_a, n, -1);
    for (int i = 0; i < n; ++i) c[i] = (ll) floor(f_a[i].x + 0.5);
  }

  void multiply (int a[], int b[], int n_a, int n_b, int c[], int mod = (int) 1e9 + 7, int dup = 0) {
    int n = n_a + n_b - 1;
    while (n != (n & -n)) n += n & -n;
    for (int i = 0; i < n; ++i) f_a[i] = f_b[i] = cplx();
    static const int magic = 15;
    for (int i = 0; i < n_a; ++i) f_a[i] = cplx(a[i] >> magic, a[i] & (1 << magic) - 1);
    for (int i = 0; i < n_b; ++i) f_b[i] = cplx(b[i] >> magic, b[i] & (1 << magic) - 1);
    dft(f_a, n, 1);
    if (dup) for (int i = 0; i < n; ++i) f_b[i] = f_a[i];
    else dft(f_b, n, 1);
    for (int i = 0; i < n; ++i) {
      int j = (n - i) % n;
      cplx x = f_a[i] + !f_a[j];
      cplx y = f_b[i] + !f_b[j];
      cplx z = !f_a[j] - f_a[i];
      cplx t = !f_b[j] - f_b[i];
      f_c[i] = (x * t + y * z) * cplx(0, 0.25);
      f_d[i] = x * y * cplx(0, 0.25) + z * t * cplx(-0.25, 0);
    }
    dft(f_c, n, -1), dft(f_d, n, -1);
    for (int i = 0; i < n; ++i) {
      ll u = ((ll) floor(f_c[i].x + 0.5)) % mod;
      ll v = ((ll) floor(f_d[i].x + 0.5)) % mod;
      ll w = ((ll) floor(f_d[i].y + 0.5)) % mod;
      c[i] = ((u << 15) + v + (w << 30)) % mod;
    }
  }

  vector <int> multiply (vector <int> a, vector <int> b, int mod = (int) 1e9 + 7) {
    static int f_a[T], f_b[T], f_c[T];
    int n_a = a.size(), n_b = b.size();
    for (int i = 0; i < n_a; ++i) f_a[i] = a[i];
    for (int i = 0; i < n_b; ++i) f_b[i] = b[i];
    multiply(f_a, f_b, n_a, n_b, f_c, mod, a == b);
    int k = n_a + n_b - 1;
    vector <int> res(k);
    for (int i = 0; i < k; ++i) res[i] = f_c[i];
    return res;
  }

  int bigMod (int a, int k, int p) {
    if (!k) return 1;
    int res = a, t = a; --k;
    while (k) {
      if (k & 1) res = (ll) res * t % p;
      t = (ll) t * t % p; k >>= 1;
    }
    return res;
  }

  vector <int> invert (vector <int> a, int n, int mod) {
    assert(a[0]);
    vector <int> x(1, bigMod(a[0], mod - 2, mod));
    while (x.size() < n) {
      vector <int> tmp(a.begin(), a.begin() + min(a.size(), 2 * x.size()));
      vector <int> nx = multiply(multiply(x, x, mod), tmp, mod);
      x.resize(2 * x.size());
      for (int i = 0; i < x.size(); ++i) {
        x[i] += x[i]; x[i] -= nx[i];
        if (x[i] < 0) x[i] += mod;
        if (x[i] >= mod) x[i] -= mod;
      }
    }
    x.resize(n); return x;
  }

  pair <vector <int>, vector <int>> divMod (vector <int> a, vector <int> b, int mod) {
    int n = a.size(), m = b.size();
    if (n < m) return make_pair(vector <int> (), a);
    reverse(a.begin(), a.end());
    reverse(b.begin(), b.end());
    vector <int> rb = invert(b, n - m + 1, mod);
    vector <int> d = multiply(a, rb, mod);
    reverse(a.begin(), a.end());
    reverse(b.begin(), b.end());
    while (d.size() > n - m + 1) d.pop_back();
    reverse(d.begin(), d.end());
    vector <int> r = multiply(d, b, mod);
    while (r.size() >= m) r.pop_back();
    for (int i = 0; i < m; ++i) {
      r[i] = a[i] - r[i];
      if (r[i] < 0) r[i] += mod;
    }
    return make_pair(d, r);
  }

  vector <int> chirp_Z_transform (vector <int> a, int z, int k, int mod) {
    int n = a.size();
    vector <int> x, y;
    int iz = bigMod(z, mod - 2, mod);
    for (int i = 0; i < n; ++i) {
      x.push_back((ll) a[i] * bigMod(z, (ll) i * i, mod) % mod);
    }
    for (int i = 1 - n; i < k; ++i) {
      y.push_back(bigMod(iz, (ll) i * i, mod));
    }
    vector <int> r = FFT::multiply(x, y, mod), res(k);
    for (int i = 0; i < k; ++i) {
      res[i] = (ll) r[i + n - 1] * bigMod(z, (ll) i * i, mod) % mod;
    }
    return res;
  }
}

const int N = 1e5 + 5;
const int MOD = (int) 1e9 + 7;

int n, m, a[N];
vector <int> st[N << 2];

void build (int p, int L, int R) {
  if (L == R) {
    st[p].resize(2);
    st[p][0] = MOD - a[L + R >> 1];
    st[p][1] = 1;
    return;
  }
  build(p << 1, L, L + R >> 1);
  build(p << 1 | 1, (L + R >> 1) + 1, R);
  st[p] = FFT::multiply(st[p << 1], st[p << 1 | 1]);
  while (st[p].size() > R - L + 2) st[p].pop_back();
}

void divide (int p, int L, int R, vector <int> poly, vector <int> &vals) {
  poly = FFT::divMod(poly, st[p], MOD).second;
  if (L == R) {
    vals[L + R >> 1] = poly[0];
    return;
  }
  divide(p << 1, L, L + R >> 1, poly, vals);
  divide(p << 1 | 1, (L + R >> 1) + 1, R, poly, vals);
}

void testDivMod() {
  int n = 1000 + rand() % 100, m = 100 + rand() % 10;
  vector <int> a, b;
  for (int i = 0; i < n; ++i) {
    a.push_back(rand() % MOD);
  }
  for (int i = 0; i < m; ++i) {
    b.push_back(rand() % MOD);
  }
  pair <vector <int>, vector <int> > res = FFT::divMod(a, b, MOD);
  vector <int> d = res.first;
  vector <int> r = res.second;
  assert(r.size() < b.size());
  vector <int> c = FFT::multiply(b, d, MOD);
  for (int i = 0; i < r.size(); ++i) {
    c[i] += r[i];
    if (c[i] >= MOD) c[i] -= MOD;
  }
  for (int i = 0; i < n; ++i) {
    assert(a[i] == c[i]);
  }
}

void testMulEval() {
  n = 5000, m = 5000;
  vector <int> poly;
  for (int i = 0; i < n; ++i) {
    a[i] = rand() % MOD;
  }
  for (int i = 0; i < m; ++i) {
    poly.push_back(rand() % MOD);
  }
  build(1, 0, n - 1);
  vector <int> vals(n);
  divide(1, 0, n - 1, poly, vals);
  for (int i = 0; i < n; ++i) {
    int t = 0;
    for (int j = (int) poly.size() - 1; j >= 0; --j) {
      t = ((ll) t * a[i] + poly[j]) % MOD;
    }
    assert(t == vals[i]);
  }
}

void test_chirp_Z_transform() {
  int n = 1234;
  vector <int> a(n);
  for (int i = 0; i < n; ++i) {
    a[i] = rand() % MOD;
  }
  int z = rand() % MOD;
  vector <int> r = FFT::chirp_Z_transform(a, z, n, MOD);
  for (int i = 0; i < n; ++i) {
    int x = FFT::bigMod(z, i + i, MOD), y = 0;
    for (int j = 0; j < n; ++j) {
      y += (ll) a[j] * FFT::bigMod(x, j, MOD) % MOD;
      if (y >= MOD) y -= MOD;
    }
    assert(y == r[i]);
  }
}

int main() {
  srand(time(0));
  testDivMod();
  testMulEval();
  test_chirp_Z_transform();
  cerr << "Correct!\n";
  cerr << "\nTime elapsed: " << 1000 * clock() / CLOCKS_PER_SEC << "ms.\n";
  return 0;
}

