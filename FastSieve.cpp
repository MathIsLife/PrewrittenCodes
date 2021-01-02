#include <bits/stdc++.h>

using namespace std;

// primes up to 5e8 within 0.35 seconds
// primes up to 1e9 within 1 second
vector <int> fastSieve (const int N, const int Q = 17, const int L = 1 << 15) {
  const int M = (N + 29) / 30;
  const int two = sqrt(N), four = sqrt(two);
  static const int r[] = {1, 7, 11, 13, 17, 19, 23, 29};
  
  struct P {
    P (int p) : p(p) {}
    int p, pos[8];
  };

  auto approxPrimeCount = [] (const int N) -> int {
    return N > 60184 ? N / (log(N) - 1.1) : max(1.0, N / (log(N) - 1.11)) + 1;
  };

  vector <bool> isPrime(two + 1, true);
  for (int i = 2; i <= four; ++i) if (isPrime[i]) {
    for (int j = i * i; j <= two; j += i) isPrime[j] = false;
  }

  const int r_size = approxPrimeCount(N + 30);

  int p_size = 3;
  vector <P> s_primes;
  vector <int> primes = {2, 3, 5};
  int p_beg = 0, prod = 1;
  primes.resize(r_size);
  
  for (int p = 7; p <= two; ++p) {
    if (!isPrime[p]) continue;
    if (p <= Q) prod *= p, ++p_beg, primes[p_size++] = p;
    auto cur = P(p);
    for (int t = 0; t < 8; ++t) {
      int j = (p <= Q) ? p : p * p;
      while (j % 30 != r[t]) j += p << 1;
      cur.pos[t] = j / 30;
    }
    s_primes.push_back(cur);
  }

  vector <unsigned char> pre(prod, 0xFF);

  for (size_t it = 0; it < p_beg; ++it) {
    auto cur = s_primes[it];
    const int p = cur.p;
    for (int t = 0; t < 8; ++t) {
      const unsigned char m = ~(1 << t);
      for (int i = cur.pos[t]; i < prod; i += p) pre[i] &= m;
    }
  }

  const int block_size = (L + prod - 1) / prod * prod;

  vector <unsigned char> block(block_size);
  unsigned char *p_block = block.data();

  for (int beg = 0; beg < M; beg += block_size, p_block -= block_size) {
    int end = min(M, beg + block_size);
    for (int i = beg; i < end; i += prod) {
      copy(pre.begin(), pre.end(), p_block + i);
    }
    if (beg == 0) p_block[0] &= 0xFE;
    for (size_t it = p_beg; it < s_primes.size(); ++it) {
      auto &cur = s_primes[it];
      const int p = cur.p;
      for (int t = 0; t < 8; ++t) {
        int i = cur.pos[t];
        const unsigned char m = ~(1 << t);
        for (; i < end; i += p) p_block[i] &= m;
        cur.pos[t] = i;
      }
    }
    for (int i = beg; i < end; ++i) {
      for (int m = p_block[i]; m > 0; m &= m - 1) {
        primes[p_size++] = i * 30 + r[__builtin_ctz(m)];
      }
    }
  }

  assert(p_size <= r_size);
  while (p_size > 0 and primes[p_size - 1] > N) --p_size;
  primes.resize(p_size); return primes;
}

int main() {
  int LIM; cin >> LIM;
  auto primes = fastSieve(LIM);
  for (int x : primes) printf("%d ", x);
  puts("");
  return 0;
}

