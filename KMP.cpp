#include <bits/stdc++.h>

using namespace std;

void buildPi (string &p, vector <int> &pi) {
  pi.resize(p.size());
  int k = -2;
  for (int i = 0; i < p.size(); ++i) {
    while (k >= -1 and p[k + 1] != p[i]) k = k == -1 ? -2 : pi[k];
    pi[i] = ++k;
  }
}

void KMP (string &t, string &p) {
  vector <int> pi;
  buildPi(p, pi);
  int k = -1;
  for (int i = 0; i < t.size(); ++i) {
    while (k >= -1 and p[k + 1] != t[i]) k = k == -1 ? -2 : pi[k];
    ++k;
    if (k == p.size() - 1) {
      cout << "Matched starting at " << i - k << '\n';
      k = k == -1 ? -2 : pi[k];
    }
  }
}

int main() {
  string text = "abacaba";
  string pattern = "aba";
  KMP(text, pattern);
  return 0;
}

