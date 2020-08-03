#include <bits/stdc++.h>

using namespace std;

typedef long long ll;

struct frac {
  ll u, d;

  frac (ll _u = 0, ll _d = 1) {
    u = _u, d = _d;
    if (d) {
      ll g = __gcd(abs(u), abs(d));
      u /= g, d /= g;
      if (d < 0) u = -u, d = -d;
    } else {
      u = 1;
    }
  }

  inline void print() {
    cerr << u << "/" << d << '\n';
  }

  frac operator + (const frac &f) const {
    return frac(u * f.d + d * f.u, d * f.d);
  }

  frac operator - (const frac &f) const {
    return frac(u * f.d - d * f.u, d * f.d);
  }

  frac operator * (const frac &f) const {
    return frac(u * f.u, d * f.d);
  }

  frac operator / (const frac &f) const {
    return frac(u * f.d, d * f.u);
  }
  
  bool operator < (const frac &f) const {
    return u * f.d < d * f.u;
  }

  bool operator <= (const frac &f) const {
    return u * f.d <= d * f.u;
  }
  
  bool operator > (const frac &f) const {
    return u * f.d > d * f.u;
  }

  bool operator >= (const frac &f) const {
    return u * f.d >= d * f.u;
  }

  bool operator == (const frac &f) const {
    return u * f.d == d * f.u;
  }
};

int main() {

  return 0;
}

