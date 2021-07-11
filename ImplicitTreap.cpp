// 0-based implicit treap with subtree sum and parent support

#include <bits/stdc++.h>

using namespace std;

typedef long long ll;

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

struct node {
  int heap, tot;
  ll value, sum;
  node *l, *r, *par;

  node (ll value) : value(value), sum(value), tot(0), heap(rng()), l(nullptr), r(nullptr), par(nullptr) {}
};

inline int size (node *t) {
  return t ? t -> tot : 0;
}

inline ll sum (node *t) {
  return t ? t -> sum : 0;
}

inline void refresh (node *t) {
  if (t) {
    t -> tot = 1 + size(t -> l) + size(t -> r);
    t -> sum = sum(t -> l) + sum(t -> r) + t -> value;
    if (t -> l) t -> l -> par = t;
    if (t -> r) t -> r -> par = t;
    t -> par = nullptr;
  }
}

void split (node *t, node* &l, node* &r, int key, int add = 0) {
  if (!t) return void(l = r = nullptr);
  int cur_key = add + size(t -> l);
  if (cur_key >= key) split(t -> l, l, t -> l, key, add), r = t;
  else split(t -> r, t -> r, r, key, cur_key + 1), l = t;
  refresh(t);
}

void merge (node* &t, node *l, node *r) {
  if (!l or !r) t = l ? l : r;
  else if (l -> heap > r -> heap) merge(l -> r, l -> r, r), t = l;
  else merge(r -> l, l, r -> l), t = r;
  refresh(t);
}

inline node* getRoot (node *u) {
  while (u -> par) u = u -> par;
  return u;
}

int main() {
  int q;
  cin >> q;
  vector <node*> a(q + 1);
  for (int i = 1; i <= q; ++i) {
    int cmd, x, y;
    scanf("%d %d", &cmd, &x);
    if (cmd == 1) {
      a[i] = new node(x);
    } else if (cmd == 2) {
      scanf("%d", &y);
      node *one = getRoot(a[x]), *two = getRoot(a[y]), *root = nullptr;
      if (one != two) merge(root, one, two);
    } else if (cmd == 3) {
      scanf("%d", &y);
      node *l, *r, *root = getRoot(a[x]);
      split(root, l, r, y);
    } else {
      printf("%lld\n", sum(getRoot(a[x])));
    }
  }
  return 0;
}


