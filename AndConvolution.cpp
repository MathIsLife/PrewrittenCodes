#include <bits/stdc++.h>

using namespace std;

void sos (vector <int> &dp, int x = 1) {
  int sz = 31 - __builtin_clz(dp.size());
  for (int i = 0; i < sz; ++i) for (int j = 0; j < 1 << sz; ++j) {
    if (j & 1 << i) dp[j ^ 1 << i] += x * dp[j];
  }
}

vector <int> andConv (vector <int> a, vector <int> b) {
  sos(a), sos(b);
  for (int i = 0; i < a.size(); ++i) a[i] *= b[i];
  sos(a, -1); return a;
}

int main() {
  
  return 0;
}

