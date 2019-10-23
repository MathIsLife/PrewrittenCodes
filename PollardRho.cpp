#include <bits/stdc++.h>

using namespace std;

namespace Rho {
  template <typename num_t>
  num_t mulmod (num_t a, num_t b, num_t p) {
    a %= p; b %= p;
    num_t q = (num_t) ((long double) a * b / p);
    num_t r = a * b - q * p;
    while (r < 0) r += p;
    while (r >= p) r -= p;
    return r;
  }

  template <typename num_t>
  num_t powmod (num_t n, num_t k, num_t p) {
    num_t r = 1;
    for (; k; k >>= 1) {
      if (k & 1) r = mulmod(r, n, p);
      n = mulmod(n, n, p);
    }
    return r;
  }

  template <typename num_t>
  int rabin (num_t n) {
    if (n == 2) return 1;
    if (n < 2 || (~n & 1)) return 0;
    const num_t p[9] = {2, 3, 5, 7, 11, 13, 17, 19, 23};
    num_t a, d = n - 1, mx = 4;
    int i, r, s = 0;
    while (~d & 1) {++s; d >>= 1;}
    for (i = 0; i < mx; ++i) {
      if (n == p[i]) return 1;
      if (!(n % p[i])) return 0;
      a = powmod(p[i], d, n);
      if (a ^ 1) {
        for (r = 0; r < s && a != n - 1; ++r) a = mulmod(a, a, n);
        if (r == s) return 0;
      }
    }
    return 1;
  }

  template <typename num_t>
  inline num_t f (num_t a, num_t b, num_t n) {
    return (mulmod(a, a, n) + b) % n;
  }

  template <typename num_t>
  void factorize (num_t n, vector <num_t> &facs) {
    static int init_seed = 0;
    if (!init_seed) {
      init_seed = 1;
      srand(2311);
    }
    if (n == 1) return;
    if (rabin(n)) return void(facs.emplace_back(n));
    if (n == 4) {
      facs.emplace_back(2);
      facs.emplace_back(2);
      return;
    }
    while (1) {
      num_t a = rand() & 63, x = 2, y = 2;
      while (1) {
        x = f(x, a, n), y = f(f(y, a, n), a, n);
        if (x == y) break;
        num_t p = __gcd(n, y <= x ? x - y : y - x);
        if (p > 1) {
          factorize(p, facs), factorize(n / p, facs);
          return;
        }
      }
    }
  }
}

int t; long long n;

int main() {
  srand(time(0));
  scanf("%lld", &n);
  vector <long long> facs;
  Rho::factorize(n, facs);
  map <long long, int> ppf;
  for (auto it : facs) ++ppf[it];
  for (auto it : ppf) cout << it.first << " " << it.second << '\n';
  return 0;
}

