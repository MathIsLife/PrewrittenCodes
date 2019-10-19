// DYNAMIC DIAMETER ONLINE
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
// ------------------------------------------------------------------------

// NUMBER THEORETIC TRANSFORM
#include <bits/stdc++.h>

using namespace std;

typedef long long ll;

const int G = 3;
const int MOD = 998244353;
const int N = (1 << 20) + 5;

int rev[N], w[N], inv_n;

int bigMod (int a, int e, int mod) {
  if (e == -1) e = mod - 2;
  int ret = 1;
  while (e) {
    if (e & 1) ret = (ll) ret * a % mod;
    a = (ll) a * a % mod; e >>= 1;
  }
  return ret;
}

void prepare (int &n) { 
  int sz = abs(31 - __builtin_clz(n));
  int r = bigMod(G, (MOD - 1) / n, MOD); 
  inv_n = bigMod(n, MOD - 2, MOD), w[0] = w[n] = 1; 
  for (int i = 1; i < n; ++i) w[i] = (ll) w[i - 1] * r % MOD;
  for (int i = 1; i < n; ++i) rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (sz - 1));
}

void ntt (int *a, int n, int dir) { 
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

int f_a[N], f_b[N];

vector <int> multiply (vector <int> a, vector <int> b) {
  int sz = 1, n = a.size(), m = b.size();
  while (sz < n + m - 1) sz <<= 1; prepare(sz);
  for (int i = 0; i < sz; ++i) f_a[i] = i < n ? a[i] : 0;
  for (int i = 0; i < sz; ++i) f_b[i] = i < m ? b[i] : 0;
  ntt(f_a, sz, 0); ntt(f_b, sz, 0);
  for (int i = 0; i < sz; ++i) f_a[i] = (ll) f_a[i] * f_b[i] % MOD;
  ntt(f_a, sz, 1); return vector <int> (f_a, f_a + n + m - 1);
}

// G = primitive_root(MOD)
int primitive_root (int p) {
  vector <int> factor;
  int tmp = p - 1;
  for (int i = 2; i * i <= tmp; ++i) {
    if (tmp % i == 0) {
      factor.emplace_back(i);
      while (tmp % i == 0) tmp /= i;  
    }
  }
  if (tmp != 1) factor.emplace_back(tmp);
  for (int root = 1; ; ++root) {
    bool flag = true;
    for (int i = 0; i < (int) factor.size(); ++i) {
      if (bigMod(root, (p - 1) / factor[i], p) == 1) {
        flag = false; break;
      }
    }
    if (flag) return root;
  }
}

int main() {
  // (x + 2)(x + 3) = x^2 + 5x + 6
  vector <int> a = {2, 1};
  vector <int> b = {3, 1};
  vector <int> c = multiply(a, b);
  for (int x : c) cout << x << " "; cout << endl;
  return 0;
}

// ------------------------------------------------------------------------

// FAST FOURIER TRANSFORM
#include <bits/stdc++.h>

using namespace std;

typedef long long ll;
typedef long double ld;

struct cplx {
  ld a, b;

  cplx (ld a = 0, ld b = 0) : a(a), b(b) {}
  
  const cplx operator + (const cplx &c) const {
    return cplx(a + c.a, b + c.b);
  }

  const cplx operator - (const cplx &c) const {
    return cplx(a - c.a, b - c.b);
  }

  const cplx operator * (const cplx &c) const {
    return cplx(a * c.a - b * c.b, a * c.b + b * c.a);
  }

  const cplx operator / (const ld &x) const {
    return cplx(a / x, b / x);
  }
};

const ld PI = acosl(-1);
const int N = (1 << 20) + 5;

int rev[N]; cplx w[N];

void prepare (int &n) { 
  int sz = abs(31 - __builtin_clz(n));
  cplx r = cplx(cosl(2 * PI / n), sinl(2 * PI / n)); 
  w[0] = w[n] = 1; 
  for (int i = 1; i < n; ++i) w[i] = w[i - 1] * r;
  for (int i = 1; i < n; ++i) rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (sz - 1));
}

void fft (cplx *a, int n, int dir) { 
  for (int i = 1; i < n - 1; ++i) { 
    if (i < rev[i]) swap(a[i], a[rev[i]]);
  }
  for (int m = 2; m <= n; m <<= 1) {
    for (int i = 0; i < n; i += m) {
      for (int j = 0; j < (m >> 1); ++j) {
        cplx &u = a[i + j], &v = a[i + j + (m >> 1)]; 
        cplx t = v * w[dir ? n - n / m * j : n / m * j];
        v = u - t, u = u + t;
      }
    }
  } if (dir) for (int i = 0; i < n; ++i) a[i] = a[i] / n;
}

cplx f_a[N], f_b[N];

vector <ll> multiply (vector <ll> a, vector <ll> b) {
  int sz = 1, n = a.size(), m = b.size();
  while (sz < n + m - 1) sz <<= 1; prepare(sz);
  for (int i = 0; i < sz; ++i) f_a[i] = i < n ? cplx(a[i]) : cplx();
  for (int i = 0; i < sz; ++i) f_b[i] = i < m ? cplx(b[i]) : cplx();
  fft(f_a, sz, 0); fft(f_b, sz, 0);
  for (int i = 0; i < sz; ++i) f_a[i] = f_a[i] * f_b[i];
  fft(f_a, sz, 1); vector <ll> ret(n + m - 1);
  for (int i = 0; i < n + m - 1; ++i) ret[i] = f_a[i].a + 0.3;
  return ret;
}

int main() {
  // (x + 2)(x + 3) = x^2 + 5x + 6
  vector <ll> a = {2, 1};
  vector <ll> b = {3, 1};
  vector <ll> c = multiply(a, b);
  for (int x : c) cout << x << " "; cout << endl;
  return 0;
}

// ------------------------------------------------------------------

// FFT APPLICATIONS
 
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
      x.emplace_back((ll) a[i] * bigMod(z, (ll) i * i, mod) % mod);
    }
    for (int i = 1 - n; i < k; ++i) {
      y.emplace_back(bigMod(iz, (ll) i * i, mod));
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
    a.emplace_back(rand() % MOD);
  }
  for (int i = 0; i < m; ++i) {
    b.emplace_back(rand() % MOD);
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
    poly.emplace_back(rand() % MOD);
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

// -------------------------------------------------------

// WALSH HADAMARD TRANSFORM
#include <bits/stdc++.h>

using namespace std;

typedef long long ll;

const int OR = 0;
const int AND = 1;
const int XOR = 2;
const int N = (1 << 20) + 5;

namespace FWHT {
  ll a[N], b[N];

  void forward_fwht (ll *arr, int n, int flag = XOR) {
    if (n == 0) return;
    int i, m = n >> 1;
    forward_fwht(arr, m, flag);
    forward_fwht(arr + m, m, flag);

    // apply mod if required
    for (i = 0; i < m; ++i) {
      ll x = arr[i], y = arr[i + m];
      if (flag == OR) arr[i] = x, arr[i + m] = x + y;
      if (flag == AND) arr[i] = x + y, arr[i + m] = y;
      if (flag == XOR) arr[i] = x + y, arr[i + m] = x - y;
    }
  }

  void inverse_fwht (ll *arr, int n, int flag = XOR) {
    if (n == 0) return;
    int i, m = n >> 1;
    inverse_fwht(arr, m, flag);
    inverse_fwht(arr + m, m, flag);

    // apply mod if required
    for (i = 0; i < m; ++i) { 
      ll x = arr[i], y = arr[i + m];
      if (flag == OR) arr[i] = x, arr[i + m] = y - x;
      if (flag == AND) arr[i] = x - y, arr[i + m] = y;
      if (flag == XOR) arr[i] = (x + y) >> 1, arr[i + m] = (x - y) >> 1;
    }
  }

  vector <ll> convolution (int n, ll *A, ll *B, int flag = XOR) {
    assert(!(n & (n - 1)));
    for (int i = 0; i < n; ++i) a[i] = A[i];
    for (int i = 0; i < n; ++i) b[i] = B[i];
    forward_fwht(a, n, flag);
    forward_fwht(b, n, flag);
    // apply mod if required
    for (int i = 0; i < n; ++i) a[i] = a[i] * b[i];
    inverse_fwht(a, n, flag);
    return vector <ll> (a, a + n);
  }
}

int n; ll A[N], B[N];

int main() {
  srand(time(0)); n = 1 << 3;
  for (int i = 0; i < n; ++i) A[i] = rand() & 1, B[i] = rand() & 1;
  for (int i = 0; i < n; ++i) cout << A[i] << " "; cout << '\n';
  for (int i = 0; i < n; ++i) cout << B[i] << " "; cout << '\n';
  vector <ll> res = FWHT::convolution(n, A, B, XOR);
  for (auto it : res) cout << it << " "; cout << '\n';
  return 0;
}

// ----------------------------------------------------

// AP FLOOR SUM

#include <bits/stdc++.h>

using namespace std;

typedef long long ll;

// sum [(x + kn) / m] for 0 <= k < m
inline ll aux (ll x, ll n, ll m) {
  ll g = __gcd(n, m);
  return g * (x / g) + (m * n - m - n + g) / 2;
}

// sum [(x + kn) / m] for 0 <= k < lim < m
ll get (ll x, ll n, ll m, ll lim) {
  if (!lim) return 0;
  ll ret = lim * (x / m) + lim * (lim - 1) * (n / m) / 2; 
  n %= m, x %= m; if (!n) return ret;
  ll p = (x + (lim - 1) * n) / m;
  return ret + p * lim - get(m - x + n - 1, m, n, p);
}

// sum [(x + kn) / m] for 0 <= k < lim in O(lg max(n, m))
// m > 0, lim >= 0
inline ll floorAPsum (ll x, ll n, ll m, ll lim) {
  ll ret = 0;
  if (x < 0) {
    ll q = x / m; x %= m;
    if (x) x += m, --q;
    ret += q * lim;
  } if (n < 0) {
    ll q = n / m; n %= m;
    if (n) n += m, --q;
    ret += lim * (lim - 1) * q / 2;
  }
  ll tmp = aux(x, n, m), tot = lim / m;
  ret += tmp * tot + tot * (tot - 1) * n * m / 2 + tot * n * (lim - tot * m);
  return ret + get(x, n, m, lim % m);
}

ll brute (ll x, ll n, ll m, ll lim) {
  ll ret = 0;
  for (ll k = 0; k < lim; ++k) {
    if ((x + k * n) % m == 0) ret += (x + k * n) / m;
    else if (x + k * n < 0) ret += (x + k * n) / m - 1;
    else ret += (x + k * n) / m;
  }
  return ret;
}

int main() {
  ll x, n, m, lim;
  while (cin >> x >> n >> m >> lim) {
    ll p = brute(x, n, m, lim);
    ll q = floorAPsum(x, n, m, lim);
    cout << p << " " << q << endl;
    assert(p == q);
  }
  return 0;
}

// -----------------------------------------------------------

// GRID PATTERN HAMMING DISTANCE

#include <bits/stdc++.h>

using namespace std;

// FFT BEGIN

typedef long long ll;
typedef long double ld;

struct cplx {
  ld a, b;

  cplx (ld a = 0, ld b = 0) : a(a), b(b) {}
  
  const cplx operator + (const cplx &c) const {
    return cplx(a + c.a, b + c.b);
  }

  const cplx operator - (const cplx &c) const {
    return cplx(a - c.a, b - c.b);
  }

  const cplx operator * (const cplx &c) const {
    return cplx(a * c.a - b * c.b, a * c.b + b * c.a);
  }

  const cplx operator / (const ld &x) const {
    return cplx(a / x, b / x);
  }
};

const ld PI = acosl(-1);
const int N = (1 << 20) + 5;

int rev[N]; cplx w[N];

void prepare (int &n) { 
  int sz = abs(31 - __builtin_clz(n));
  cplx r = cplx(cosl(2 * PI / n), sinl(2 * PI / n)); 
  w[0] = w[n] = 1; 
  for (int i = 1; i < n; ++i) w[i] = w[i - 1] * r;
  for (int i = 1; i < n; ++i) rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (sz - 1));
}

void fft (cplx *a, int n, int dir) { 
  for (int i = 1; i < n - 1; ++i) { 
    if (i < rev[i]) swap(a[i], a[rev[i]]);
  }
  for (int m = 2; m <= n; m <<= 1) {
    for (int i = 0; i < n; i += m) {
      for (int j = 0; j < (m >> 1); ++j) {
        cplx &u = a[i + j], &v = a[i + j + (m >> 1)]; 
        cplx t = v * w[dir ? n - n / m * j : n / m * j];
        v = u - t, u = u + t;
      }
    }
  } if (dir) for (int i = 0; i < n; ++i) a[i] = a[i] / n;
}

cplx f_a[N], f_b[N];

vector <ll> multiply (vector <ll> a, vector <ll> b) {
  int sz = 1, n = a.size(), m = b.size();
  while (sz < n + m - 1) sz <<= 1; prepare(sz);
  for (int i = 0; i < sz; ++i) f_a[i] = i < n ? cplx(a[i]) : cplx();
  for (int i = 0; i < sz; ++i) f_b[i] = i < m ? cplx(b[i]) : cplx();
  fft(f_a, sz, 0); fft(f_b, sz, 0);
  for (int i = 0; i < sz; ++i) f_a[i] = f_a[i] * f_b[i];
  fft(f_a, sz, 1); vector <ll> ret(n + m - 1);
  for (int i = 0; i < n + m - 1; ++i) ret[i] = f_a[i].a + 0.3;
  return ret;
}

// FFT END

const int M = 1010;

char g[M][M], pat[M][M];
int n, m, r, c, res[M][M];

vector <int> strVec, patternVec, mulVec;

int match() {
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      res[i][j] = r * c;
    }
  }
  for (int ch = 'a'; ch <= 'z'; ++ch) {
    int cnt = 0;
    vector <ll> strVec, patVec;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < m; ++j) {
        if (i < r and j < c and pat[i][j] == ch) {
          ++cnt;
          patVec.emplace_back(1);
        } else {
          patVec.emplace_back(0);
        }
        strVec.emplace_back(g[i][j] == ch);
      }
    }
    reverse(patVec.begin(), patVec.end());
    vector <ll> prod = multiply(strVec, patVec);
    for (int i = 0; i <= n - r; ++i) {
      for (int j = 0; j <= m - c; ++j) {
        if (i * m + j + m * n - 1 < prod.size()) {
          res[i][j] -= prod[i * m + j + m * n - 1];
        }
      }
    }
  }
  // res[i][j] stores the hamming distance
  int ret = 0;
  for (int i = 0; i <= n - r; ++i) {
    for (int j = 0; j <= m - c; ++j) ret += res[i][j] == 0;
  }
  return ret;
}

int main() {
  cin >> n >> m >> r >> c;
  for (int i = 0; i < n; ++i) scanf("%s", g[i]);
  for (int i = 0; i < r; ++i) scanf("%s", pat[i]);
  cout << match() << '\n';
  return 0;
}

// -------------------------------------------------

// MAXIMUM PARTITION SORTED SUM
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

// --------------------------------------------------

// MINIMUM DIAMETER SPANNING TREE
#include <bits/stdc++.h>

using namespace std;

typedef long long ll;
typedef __int128 lll;

const int N = 505;
const int M = 200010;
const ll INF = 1e18 + 5;

inline bool good (lll x_1, lll y_1, lll x_2, lll y_2, lll x_3, lll y_3) {
  return x_1 * y_2 + x_2 * y_3 + x_3 * y_1 >= y_1 * x_2 + y_2 * x_3 + y_3 * x_1;
}

ll d[N][N], W[M];
vector <int> g[N];
int n, m, U[M], V[M], id[N], par[N];
ll ed, far, global = LLONG_MAX, dist[N];

int main() {
  cin >> n >> m;
  memset(d, 63, sizeof d);
  for (int i = 1; i <= n; ++i) {
    d[i][i] = 0, id[i] = par[i] = i;
  }
  for (int i = 1; i <= m; ++i) {
    scanf("%d %d %lld", U + i, V + i, W + i);
    if (U[i] > V[i]) swap(U[i], V[i]);
    W[i] <<= 1LL;
    int u = U[i], v = V[i]; ll w = W[i];
    d[u][v] = min(d[u][v], w);
    d[v][u] = min(d[v][u], w);
    g[u].emplace_back(i), g[v].emplace_back(i);
  }
  for (int k = 1; k <= n; ++k) {
    for (int i = 1; i <= n; ++i) {
      for (int j = 1; j <= n; ++j) {
        d[i][j] = min(d[i][j], d[i][k] + d[k][j]);
      }
    }
  }
  for (int u = 1; u <= n; ++u) {
    sort(id + 1, id + n + 1, [&] (int i, int j) {return d[u][i] > d[u][j];});
    for (int e : g[u]) {
      int v = U[e] ^ u ^ V[e]; ll w = W[e];
      if (v < u) continue; 
      int last = id[1];
      ll opt = 0, minVal = d[u][last];
      for (int it = 2; it <= n; ++it) {
        int i = id[it];
        if (d[v][i] <= d[v][last]) continue; 
        ll curX = (w - d[u][i] + d[v][last]) / 2;
        ll curY = (w + d[u][i] + d[v][last]) / 2;
        if (minVal > curY) minVal = curY, opt = curX; 
        last = i;
      }
      if (minVal > d[v][last]) minVal = d[v][last], opt = w;
      if (minVal < global) global = minVal, ed = e, far = opt;
    }
  }
  ++n;
  ++m, U[m] = U[ed], V[m] = n, W[m] = far;
  g[U[m]].emplace_back(m), g[V[m]].emplace_back(m);
  ++m, U[m] = V[ed], V[m] = n, W[m] = W[ed] - far; 
  g[U[m]].emplace_back(m), g[V[m]].emplace_back(m);
  priority_queue <pair <ll, int>> pq;
  pq.emplace(0, n);
  for (int i = 1; i < n; ++i) dist[i] = INF;
  while (!pq.empty()) {
    int u = pq.top().second; pq.pop();
    for (int e : g[u]) if (e ^ ed) {
      int v = U[e] ^ u ^ V[e]; ll w = W[e];
      if (dist[v] > dist[u] + w) {
        dist[v] = dist[u] + w, par[v] = u, pq.emplace(-dist[v], v);
      }
    }
  }
  printf("%lld\n", global);
  printf("%d %d\n", U[ed], V[ed]);
  for (int i = 1; i < n; ++i) if (par[i] ^ n) {
    int u = i, v = par[i];
    if (u > v) swap(u, v); 
    printf("%d %d\n", u, v);
  }
  return 0;
}

// -------------------------------------------------------

// PUSH RELABEL
// O(V^2 sqrt(E)), solves SPOJ FASTFLOW

#include <bits/stdc++.h>

using namespace std;

typedef long long ll;

struct edge {
  int dest, back; ll f, c;

  edge (int d, int b, ll f, ll c) : dest(d), back(b), f(f), c(c) {}
};

struct PushRelabel {
  vector <vector <edge>> g;
  vector <ll> ec;
  vector <edge*> cur;
  vector <vector <int>> hs; 
  vector <int> H;
  
  PushRelabel (int n) : g(n), ec(n), cur(n), hs(n + n), H(n) {}

  void AddEdge (int s, int t, ll cap, ll rcap = 0) {
    if (s == t) return;
    g[s].emplace_back(t, g[t].size(), 0, cap);
    g[t].emplace_back(s, g[s].size() - 1, 0, rcap);
  }

  void AddFlow (edge &e, ll f) {
    edge &back = g[e.dest][e.back];
    if (!ec[e.dest] and f) hs[H[e.dest]].emplace_back(e.dest);
    e.f += f; e.c -= f; ec[e.dest] += f;
    back.f -= f; back.c += f; ec[back.dest] -= f;
  }

  ll MaxFlow (int s, int t) {
    int v = g.size(); H[s] = v; ec[t] = 1;
    vector <int> co(v + v); co[0] = v - 1;
    for (int i = 0; i < v; ++i) cur[i] = g[i].data();
    for (auto &e : g[s]) AddFlow(e, e.c);
    for (int hi = 0; ;) {
      while (hs[hi].empty()) if (!hi--) return -ec[s];
      int u = hs[hi].back(); hs[hi].pop_back();
      while (ec[u] > 0) 
        if (cur[u] == g[u].data() + g[u].size()) {
          H[u] = 1e9 + 5;
          for (auto &e : g[u]) if (e.c and H[u] > H[e.dest] + 1)
            H[u] = H[e.dest] + 1, cur[u] = &e;
          if (++co[H[u]], !--co[hi] and hi < v)
            for (int i = 0; i < v; ++i) if (hi < H[i] and H[i] < v)
              --co[H[i]], H[i] = v + 1;
          hi = H[u];
        } else if (cur[u] -> c and H[u] == H[cur[u] -> dest] + 1)
          AddFlow(*cur[u], min(ec[u], cur[u] -> c));
        else ++cur[u];
    }
  }
};

int main() {
  int n, m;
  cin >> n >> m;
  PushRelabel flow(n);
  while (m--) {
    int u, v, c;
    scanf("%d %d %d", &u, &v, &c);
    --u, --v;
    flow.AddEdge(u, v, c, c);
  }
  cout << flow.MaxFlow(0, n - 1) << '\n';
  return 0;
}

// ---------------------------------------------------------

// MINIMUM ENCLOSING CIRCLE
// Expected runtime: O(n)
// Solves Gym 102299J

#include <bits/stdc++.h>

using namespace std;

typedef long double ld;
typedef pair <ld, ld> point;

#define x first
#define y second

point operator + (const point &a, const point &b) {
  return point(a.x + b.x, a.y + b.y);
}

point operator - (const point &a, const point &b) {
  return point(a.x - b.x, a.y - b.y);
}

point operator * (const point &a, const ld &b) {
  return point(a.x * b, a.y * b);
}

point operator / (const point &a, const ld &b) {
  return point(a.x / b, a.y / b);
}

const ld EPS = 1e-8;
const ld INF = 1e20;
const ld PI = acosl(-1);

inline ld dist (point a, point b) {
  return hypotl(a.x - b.x, a.y - b.y);
}

inline ld sqDist (point a, point b) {
  return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
}

inline ld dot (point a, point b) {
  return a.x * b.x + a.y * b.y;
}

inline ld cross (point a, point b) {
  return a.x * b.y - a.y * b.x;
}

inline ld cross (point a, point b, point c) {
  return cross(b - a, c - a);
}

inline point perp (point a) {
  return point(-a.y, a.x);
}

// circle through 3 points
pair <point, ld> getCircle (point a, point b, point c) {
  pair <point, ld> ret;
  ld den = (ld) 2 * cross(a, b, c);
  ret.x.x = ((c.y - a.y) * (dot(b, b) - dot(a, a)) - (b.y - a.y) * (dot(c, c) - dot(a, a))) / den;
  ret.x.y = ((b.x - a.x) * (dot(c, c) - dot(a, a)) - (c.x - a.x) * (dot(b, b) - dot(a, a))) / den;
  ret.y = dist(ret.x, a);
  return ret;
}

pair <point, ld> minCircleAux (vector <point> &s, point a, point b, int n) {
  ld lo = -INF, hi = INF;
  for (int i = 0; i < n; ++i) {
    auto si = cross(b - a, s[i] - a);
    if (fabs(si) < EPS) continue;
    point m = getCircle(a, b, s[i]).x;
    auto cr = cross(b - a, m - a);
    si < 0 ? hi = min(hi, cr) : lo = max(lo, cr);
  }
  ld v = 0 < lo ? lo : hi < 0 ? hi : 0;
  point c = (a + b) * 0.5 + perp(b - a) * v / sqDist(a, b);
  return {c, sqDist(a, c)};
}

pair <point, ld> minCircle (vector <point> &s, point a, int n) {
  random_shuffle(s.begin(), s.begin() + n);
  point b = s[0], c = (a + b) * 0.5;
  ld r = sqDist(a, c);
  for (int i = 1; i < n; ++i) {
    if (sqDist(s[i], c) > r * (1 + EPS)) {
      tie(c, r) = n == s.size() ? minCircle(s, s[i], i) : minCircleAux(s, a, s[i], i);
    }
  }
  return {c, r};
}

pair <point, ld> minCircle (vector <point> s) {
  assert(!s.empty()); 
  if (s.size() == 1) return {s[0], 0};
  return minCircle(s, s[0], s.size());
}

int n; vector <point> p;

int main() {
  cin >> n;
  while (n--) {
    double x, y;
    scanf("%lf %lf", &x, &y);
    p.emplace_back(x, y);
  }
  pair <point, ld> circ = minCircle(p);
  printf("%0.12f %0.12f %0.12f\n", (double) circ.x.x, (double) circ.x.y, (double) (0.5 * circ.y));
  return 0;
}

// ---------------------------------------------------





