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

  const cplx conj() const {
    return cplx(a, -b);
  }
};

const ld PI = acosl(-1);
const int N = (1 << 20) + 5;

int rev[N]; cplx w[N];

void prepare (int &n) { 
  int sz = __builtin_ctz(n);
  for (int i = 1; i < n; ++i) rev[i] = (rev[i >> 1] >> 1) | ((i & 1) << (sz - 1));
  w[0] = 0, w[1] = 1, sz = 1; 
  while (1 << sz < n) {
    cplx w_n = cplx(cosl(2 * PI / (1 << (sz + 1))), sinl(2 * PI / (1 << (sz + 1))));
    for (int i = 1 << (sz - 1); i < (1 << sz); ++i) {
      w[i << 1] = w[i], w[i << 1 | 1] = w[i] * w_n; 
    } ++sz;
  }
}

void fft (cplx *a, int n) { 
  for (int i = 1; i < n - 1; ++i) { 
    if (i < rev[i]) swap(a[i], a[rev[i]]);
  }
  for (int h = 1; h < n; h <<= 1) {
    for (int s = 0; s < n; s += h << 1) {
      for (int i = 0; i < h; ++i) {
        cplx &u = a[s + i], &v = a[s + i + h], t = v * w[h + i];
        v = u - t, u = u + t;
      }
    }
  }
}

static cplx f[N];

vector <ll> multiply (vector <ll> a, vector <ll> b) {
  int n = a.size(), m = b.size(), sz = 1;
  while (sz < n + m - 1) sz <<= 1; prepare(sz);
  for (int i = 0; i < sz; ++i) f[i] = cplx(i < n ? a[i] : 0, i < m ? b[i] : 0);
  fft(f, sz);
  for (int i = 0; i <= (sz >> 1); ++i) {
    int j = (sz - i) & (sz - 1);
    cplx x = (f[i] * f[i] - (f[j] * f[j]).conj()) * cplx(0, -0.25);
    f[j] = x, f[i] = x.conj();
  }
  fft(f, sz); vector <ll> c(n + m - 1);
  for (int i = 0; i < n + m - 1; ++i) c[i] = f[i].a / sz + 0.3;
  return c;
}

int main() {
  // (x + 2)(x + 3) = x^2 + 5x + 6
  vector <ll> a = {2, 1};
  vector <ll> b = {3, 1};
  vector <ll> c = multiply(a, b);
  for (int x : c) cout << x << " "; cout << endl;
  return 0;
}

