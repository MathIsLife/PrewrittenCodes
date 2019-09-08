#include <bits/stdc++.h>

using namespace std;

typedef long long ll;

const int OR = 0;
const int AND = 1;
const int XOR = 2;
const int N = (1 << 20) + 5;

namespace FWHT {
  ll a[N], b[N];

  void forward_fwht (ll *arr, int n, int flag = XOR) {
    if (n == 0) return;
    int i, m = n >> 1;
    forward_fwht(arr, m, flag);
    forward_fwht(arr + m, m, flag);

    // apply mod if required
    for (i = 0; i < m; ++i) {
      ll x = arr[i], y = arr[i + m];
      if (flag == OR) arr[i] = x, arr[i + m] = x + y;
      if (flag == AND) arr[i] = x + y, arr[i + m] = y;
      if (flag == XOR) arr[i] = x + y, arr[i + m] = x - y;
    }
  }

  void inverse_fwht (ll *arr, int n, int flag = XOR) {
    if (n == 0) return;
    int i, m = n >> 1;
    inverse_fwht(arr, m, flag);
    inverse_fwht(arr + m, m, flag);

    // apply mod if required
    for (i = 0; i < m; ++i) { 
      ll x = arr[i], y = arr[i + m];
      if (flag == OR) arr[i] = x, arr[i + m] = y - x;
      if (flag == AND) arr[i] = x - y, arr[i + m] = y;
      if (flag == XOR) arr[i] = (x + y) >> 1, arr[i + m] = (x - y) >> 1;
    }
  }

  vector <ll> convolution (int n, ll *A, ll *B, int flag = XOR) {
    assert(!(n & (n - 1)));
    for (int i = 0; i < n; ++i) a[i] = A[i];
    for (int i = 0; i < n; ++i) b[i] = B[i];
    forward_fwht(a, n, flag);
    forward_fwht(b, n, flag);
    for (int i = 0; i < n; ++i) a[i] = a[i] * b[i];
    inverse_fwht(a, n, flag);
    return vector <ll> (a, a + n);
  }
}

int n; ll A[N], B[N];

int main() {
  srand(time(0)); n = 1 << 3;
  for (int i = 0; i < n; ++i) A[i] = rand() & 1, B[i] = rand() & 1;
  for (int i = 0; i < n; ++i) cout << A[i] << " "; cout << '\n';
  for (int i = 0; i < n; ++i) cout << B[i] << " "; cout << '\n';
  vector <ll> res = FWHT::convolution(n, A, B, XOR);
  for (auto it : res) cout << it << " "; cout << '\n';
  return 0;
}

