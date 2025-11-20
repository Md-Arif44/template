
/**
 * Author: Arman Ferdous
 * Date:
 * License:
 * Source:
 * Description: Static hashing for 0-indexed string. Intervals are $[l, r]$.
 * Time: 
 * Status: Tested
 */

template<const ll M, const ll B> 
struct Hashing {
  int n; V<ll> h, pw;
  Hashing(const string &s) : n(sz(s)),h(n+1),pw(n+1) {
    pw[0] = 1; // ^^ s is 0 indexed
    for (int i = 1; i <= n; ++i)
      pw[i] = (pw[i-1] * B) % M,
      h[i] = (h[i-1] * B + s[i-1]) % M;
  }
  ll eval(int l, int r) { // assert(l <= r);
    return (h[r+1] - ((h[l] * pw[r-l+1])%M) + M)%M;
} };
struct Double_Hash {
  using H1 = Hashing<916969619, 101>;
  using H2 = Hashing<285646799, 103>;
  H1 h1; H2 h2;
  Double_Hash(const string &s):h1(s),h2(s){}
  pii eval(int l, int r) 
    { return {h1.eval(l,r), h2.eval(l,r)}; }
};

template <typename T> 
vector<int> z_function(const T &s) {
  int n= (int)s.size() ;
  vector<int> z(n, n);
  int l = 0, r = 0;
  for (int i = 1; i < n; i++) {
    z[i] = (i > r ? 0 : min(r - i + 1, z[i - l]));
    while (i + z[i] < n && s[z[i]] == s[i + z[i]]) {
      z[i]++;
    }
    if (i + z[i] - 1 > r) {
      l = i;
      r = i + z[i] - 1;
    }
  }
  return z;
}

template <typename T>
vector<int> Kmp(const T &s) {
  int n=(int)s.size() ;
    vector<int> p(n, 0);
      int k = 0;
      for (int i = 1; i < n; i++) {
        while (k > 0 && !(s[i] == s[k])) {    
          k = p[k - 1];
        }
        if (s[i] == s[k]) {
          k++;
        }
        p[i] = k;
      }
    return p;
}

/**
 * Author: Arman Ferdous
 * Date:
 * License:
 * Source: Arman
 * Description: Hashing with point updates on string (0-indexed). \texttt{upd(i, x): s[i] += x}. Intervals are $[l,r]$.

 * Time: $O(n \log n)$
 * Status: Tested
 */
// update delta  = new c-  s[pos ]; 
template<const ll M, const ll B> 
struct Dynamic_Hashing {
  int n; V<ll> h, pw;
  void upd(int pos, int c_add) {
    if (c_add < 0) c_add = (c_add + M) % M;
    for (int i = ++pos; i <= n; i += i&-i)
      h[i] = (h[i]+c_add *1LL* pw[i - pos])%M;
  }
  ll get(int pos, int r = 0) {
    for (int i = ++pos, j = 0; i; i -= i&-i) {
      r = (r + h[i] * 1LL * pw[j]) % M;
      j += i&-i;
    } return r;
  }
  Dynamic_Hashing(const string &s) : n(sz(s)), h(n+1), pw(n+1) {
    pw[0] = 1; // ^^ s is 0 indexed
    for (int i = 1; i <= n; ++i) pw[i] = (pw[i-1] * 1LL * B) % M;
    for (int i = 0; i < n; ++i) upd(i, s[i]);
  }
  ll eval(int l, int r) { // assert(l <= r);
    return (get(r) - ((get(l-1) * 1LL * pw[r-l+1]) % M) + M) % M;
} };
struct Double_Dynamic {
  using DH1 = Dynamic_Hashing<916969619, 571>;
  using DH2 = Dynamic_Hashing<285646799, 953>;
  DH1 h1; DH2 h2;
  Double_Dynamic(const string &s) : h1(s), h2(s) {}
  void upd(int pos, int c_add) {
    h1.upd(pos, c_add);
    h2.upd(pos, c_add);
  } pll eval(int l, int r) 
    { return {h1.eval(l,r), h2.eval(l,r)}; }
};
 
