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

struct Dsu{
   vector<int>pa,sz;
    Dsu(int n){ 
      pa.resize(n);
      sz=vector<int>(n,1);
      iota(pa.begin(),pa.end(),0);
    }
    int get(int a){ 
      return pa[a]=((pa[a]==a)?a:get(pa[a]));
    }
    bool add(int a,int b){
       a=get(a),b=get(b);
       if(a!=b){
         if(sz[a]<sz[b])swap(a,b);
          
           sz[a]+=sz[b]; pa[b]=a; 
        return true; 
       }
       else return false;
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


struct Trie{
   const int B=30;
   struct node{
       node *nxt[2] ;
       int sz;
       node(){nxt[0]=nxt[1]=NULL ;sz=0 ;}
   }*root;
   
   Trie (){ root=new node(); }
   void insert(int x){
     node *cur= root ;
      cur->sz++ ;
      for(int i=B;i>=0;i--){
          int b=x>>i&1;
          if(cur->nxt[b]==NULL)cur->nxt[b]=new node() ,cur=cur->nxt[b];
          else cur=cur->nxt[b] ;
          cur->sz++;
      }
   }
   ll get_max(int x ){
      node *cur= root ;
       ll ans =0 ;
      for(int i=B;i>=0;i--){
          int b=x>>i&1;
          if(cur->nxt[!b] && cur->nxt[!b]->sz>0)cur=cur->nxt[!b] ,ans+=(1<<i);
          else cur=cur->nxt[b] ;
      }
      return ans ; 
   }
   int query(int x, int k) { // number of values s.t. val ^ x < k
    node* cur = root;
     int  ans = 0;
    for (int i = B ; i >= 0; i--) {
      if (cur == NULL) break;
      int b1 = x >> i & 1, b2 = k >> i & 1;
      if (b2 == 1) {
        if (cur -> nxt[b1]) ans += cur -> nxt[b1] -> sz;
        cur = cur -> nxt[!b1];
      } else cur = cur -> nxt[b1];
    }
    return ans;
  }
  void remove(int val) {
    node* cur = root;
    cur -> sz--;
    for (int i = B; i >= 0; i--) {
        int b = val >> i & 1;
        cur = cur -> nxt[b];
        cur -> sz--;
    }
  }
  void del(node* cur) {
    for (int i=0;i<2;i++) if(cur->nxt[i])del(cur->nxt[i]);
    delete(cur);
  }  
};

using TY=  int  ;
using Mat=vector<vector<TY>> ;
 
Mat operator*(const Mat& a, const Mat& b) {
  if (a.empty() || b.empty())return {{}};
  
  Mat c(sz(a) , V<TY>(sz(b[0])) );
 
  rep(i,sz(c)){
    rep(j,sz(c[0])) {
      c[i][j] = 0;
      rep(k,sz(b)) {
        c[i][j] += a[i][k] * b[k][j];
      }
    }
  }
  return c;
}
 
Mat &operator*=(Mat &a, const Mat& b) {
  return a = a * b;
}
 
Mat power(Mat a,  ll b) {
  assert(b >= 0);
  Mat res(sz(a),V<TY>(sz(a[0])) );
 
  rep(i,sz(a))res[i][i] = 1;
   
  for(;b;a*=a,b>>=1)if(b&1)res*=a;
  
  return res;
}

#include<bits/stdc++.h>
using namespace std;
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
template<class T>using ordered_set=tree<T, null_type, less<T>, rb_tree_tag,tree_order_statistics_node_update>;
template<class T>using ordered_multiset=tree<T, null_type,less_equal<T>, rb_tree_tag,tree_order_statistics_node_update>;

/**
 * PBDS Tips:
 * - less<T>         : Sorted set (asc)
 * - less_equal<T>   : Multiset (asc)
 * - greater<T>      : Sorted set (desc)
 * - greater_equal<T>: Multiset (desc)
 *
 * name.order_of_key(k)  -> Count of elements < k
 * *name.find_by_order(k)-> k-th smallest element (0-indexed)
*/
