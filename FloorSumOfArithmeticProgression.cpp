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

