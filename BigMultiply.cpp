#include <bits/stdc++.h>

using namespace std;

typedef long long ll;
typedef long double ld;
typedef unsigned long long ull;

// a and b can be long long 
inline ull bigMultiply (ull a, ull b, ull M) {
  ll ret = a * b - M * ull(ld(a) * ld(b) / ld(M));
  return ret < 0 ? ret + M : (ret >= ll(M) ? ret - M : ret);
}

int main() {
  
  return 0;
}

