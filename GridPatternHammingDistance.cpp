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

