#include <bits/stdc++.h>

using namespace std;

typedef long long ll;

inline ll sum (ll n) {
  return n * (n - 1) / 2;
}

// sum [(a + di) / m] for 0 <= i < n
ll floorAPsum (ll a, ll d, ll m, ll n) {
  ll res = d / m * sum(n) + a / m * n;
  d %= m, a %= m; if (!d) return res;
  ll to = (n * d + a) / m;
  return res + (n - 1) * to - floorAPsum(m - 1 - a, m, d, to);
}

// sum (a + di) % m for 0 <= i < n
ll modAPsum (ll a, ll d, ll m, ll n) {
  a = ((a % m) + m) % m, d = ((d % m) + m) % m;
  return n * a + d * sum(n) - m * floorAPsum(a, d, m, n);
}

int main() {
  int t; 
  cin >> t;
  while (t--) {
    int n, m, a, b;
    scanf("%d %d %d %d", &n, &m, &a, &b);
    printf("%lld\n", floorAPsum(b, a, m, n));
  }
  return 0;
}

