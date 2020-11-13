#include <bits/stdc++.h>

using namespace std;

typedef long long ll;

__attribute__((target("avx"), optimize("O3", "unroll-loops")))
ll prime_pi (const ll n) {
  if (n <= 1) return 0;
  if (n == 2) return 1;
  const int sq = sqrtl(n);
  int s = 1 + sq >> 1;
  vector <int> smalls(s); for (int i = 1; i < s; ++i) smalls[i] = i;
  vector <int> roughs(s); for (int i = 0; i < s; ++i) roughs[i] = i << 1 | 1;
  vector <ll> larges(s); for (int i = 0; i < s; ++i) larges[i] = (n / (i << 1 | 1) - 1) >> 1;
  vector <bool> skip(sq + 1);
  const auto divide = [] (ll n, ll d) -> int {return (long double) n / d;};
  const auto half = [] (int n) -> int {return (n - 1) >> 1;};
  int cnt = 0;
  for (int p = 3; p <= sq; p += 2) if (!skip[p]) {
    int q = p * p;
    if ((ll) q * q > n) break; skip[p] = 1;
    for (int i = q; i <= sq; i += p << 1) skip[i] = 1;
    int ptr = 0;
    for (int k = 0; k < s; ++k) {
      int i = roughs[k];
      if (skip[i]) continue;
      ll d = (ll) i * p;
      larges[ptr] = larges[k] - (d <= sq ? larges[smalls[d >> 1] - cnt] : smalls[half(divide(n, d))]) + cnt;
      roughs[ptr++] = i;
    }
    s = ptr;
    for (int i = half(sq), j = ((sq / p) - 1) | 1; j >= p; j -= 2) {
      int c = smalls[j >> 1] - cnt;
      for (int e = (j * p) >> 1; i >= e; --i) smalls[i] -= c;
    }
    ++cnt;
  }
  larges[0] += (ll) (s + 2 * (cnt- 1)) * (s - 1) / 2;
  for (int k = 1; k < s; ++k) larges[0] -= larges[k];
  for (int l = 1; l < s; ++l) {
    int q = roughs[l]; ll m = n / q;
    int e = smalls[half(m / q)] - cnt;
    if (e < l + 1) break; ll t = 0;
    for (int k = l + 1; k <= e; ++k) t += smalls[half(divide(m, roughs[k]))];
    larges[0] += t - (ll) (e - l) * (cnt + l - 1);
  }
  return larges[0] + 1;
}

int main() {
  ll n; 
  cin >> n;
  cout << prime_pi(n) << '\n';
  return 0;
}

