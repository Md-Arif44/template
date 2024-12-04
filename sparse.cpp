template <typename T, typename F>
struct DisjointSparseTable {
  int n;
  vector<vector<T>> mat;
  F func;
  DisjointSparseTable(const vector<T>& a, const F& f) : n(int(a.size())), func(f) {
    mat.push_back(a);
    for (int p = 1; (1 << p) < n; p++) {
      mat.emplace_back(n);
      for (int mid = 1 << p; mid < n; mid += 1 << (p + 1)) {
        mat[p][mid - 1] = a[mid - 1];
        for (int j = mid - 2; j >= mid - (1 << p); j--) {
          mat[p][j] = func(a[j], mat[p][j + 1]);
        }
        mat[p][mid] = a[mid];
        for (int j = mid + 1; j < min(n, mid + (1 << p)); j++) {
          mat[p][j] = func(mat[p][j - 1], a[j]);
        }
      }
    }
  }

  T Query(int l, int r) const {
    assert(0 <= l && l <= r && r <= n-1);
    if (r - l == 0) {
      return mat[0][l];
    }
    int p = bit_width(unsigned(l^r))-1;
    return func(mat[p][l], mat[p][r]);
  }
};
template <typename T, typename F>
struct SparseTable {
  int n; //min,max,gcd
  vector<vector<T>> mat;
  F func;
  SparseTable(const vector<T>& a, const F& f) : func(f) {
    n = static_cast<int>(a.size());
    int max_log = 32 - __builtin_clz(n);
    mat.resize(max_log);
    mat[0] = a;
    for (int j = 1; j < max_log; j++) {
      mat[j].resize(n - (1 << j) + 1);
      for (int i = 0; i <= n - (1 << j); i++) {
        mat[j][i] = func(mat[j - 1][i], mat[j - 1][i + (1 << (j - 1))]);
      }
    }
  }
  T get(int from, int to) const {
    assert(0 <= from && from <= to && to <= n - 1);
    int lg = 32 - __builtin_clz(to - from + 1) - 1;
    return func(mat[lg][from], mat[lg][to - (1 << lg) + 1]);
  }
};
template <typename T, typename F>
struct Queue {
  vector<T> pref;
  vector<pair<T,T>> suf;
  F func;
  Queue(const F& f) : func(f) {}

  bool Empty() { return pref.empty() && suf.empty(); }
  int Size() { return int(pref.size()) + int(suf.size()); }
  void Clear() { pref.clear(); suf.clear(); }

  void Push(T t) {
    if (suf.empty()) {
      suf.emplace_back(t, t);
    } else {
      suf.emplace_back(t, func(suf.back().second, t));
    }
  }

  void Pop() {
    if (!pref.empty()) {
      pref.pop_back();
      return;
    }
    assert(!suf.empty());
    if (suf.size() > 1) {
      pref.resize(suf.size() - 1);
      pref[0] = suf.back().first;
      for (int i = 1; i < int(pref.size()); i++) {
        pref[i] = func(suf[int(suf.size()) - 1 - i].first, pref[i - 1]);
      }
    }
    suf.clear();
  }

  T Get() {
    assert(!Empty());
    if (pref.empty()) {
      return suf.back().second;
    }
    if (suf.empty()) {
      return pref.back();
    }
    return func(pref.back(), suf.back().second);
  
  }
};
