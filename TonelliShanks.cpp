#include <bits/stdc++.h>

using namespace std;

typedef long long ll;

int MOD;

ll bigMod (ll a, ll e) {
  if (e == -1) e = MOD - 2;
  ll ret = 1;
  while (e) {
    if (e & 1) ret = ret * a % MOD;
    a = a * a % MOD, e >>= 1;
  }
  return ret;
}

// 0 <= x < MOD
// returns sqrt(x) modulo MOD (prime), or -1 if doesn't exist
// some of the code can be made global
ll tonelliShanks (ll x) {
  if (!x) return x;
  if (MOD == 2) return x;
  if (bigMod(x, (MOD - 1) >> 1) != 1) return -1;
  ll q = MOD - 1, s = 0, z = 1;
  while (~q & 1) q >>= 1, ++s;
  if (s == 1) return bigMod(x, 1 + MOD >> 2);
  while (bigMod(z, (MOD - 1) >> 1) == 1) z = 1 + rand() % (MOD - 1);
  ll c = bigMod(z, q), r = bigMod(x, 1 + q >> 1), t = bigMod(x, q), m = s;
  while (t != 1) {
    ll i = 0, cur = t;
    while (cur != 1) cur *= cur, cur %= MOD, ++i;
    ll b = bigMod(c, 1LL << (m - i - 1));
    r *= b, r %= MOD, c = b * b % MOD, t *= c, t %= MOD, m = i;
  }
  return r;
}

int t; ll x;

int main() {
  cin >> t;
  while (t--) {
    scanf("%lld %d", &x, &MOD);
    printf("%lld\n", tonelliShanks(x));
  }
  return 0;
}

