// template 
#include<bits/stdc++.h> 
using namespace std;
 
#define TT template <typename T
 
TT,typename=void> struct cerrok:false_type {};
TT> struct cerrok <T, void_t<decltype(cerr << declval<T>() )>> : true_type {};
 
TT> constexpr void p1 (const T &x);
TT, typename V> void p1(const pair<T, V> &x) {
  cerr << "{"; p1(x.first); cerr << ","; 
  p1(x.second); cerr << "}"; 
}
TT> constexpr void p1 (const T &x) {
  if constexpr (cerrok<T>::value) cerr << x;
  else { int f = 0; cerr << '{';
    for (auto &i: x) 
      cerr << (f++ ? "," : ""), p1(i);
    cerr << "}";
} }
void p2() { cerr << "]\n"; }
TT, typename... V> void p2(T t, V... v) {
  p1(t);
  if (sizeof...(v)) cerr << ",";
  p2(v...);
}
 
#ifndef ONLINE_JUDGE 
#define dbg(x...) {cerr<<__LINE__<<"[" << #x << "]=["; p2(x); cerr << "\e ";}
#else
#define dbg(...)
#endif
 
typedef long long ll;
template<class T> using V=vector<T>;
#define all(x) (x).begin(), (x).end()
 
#define rall(x) (x).rbegin(), (x).rend()
#define ff first
#define ss second 
#define pb push_back  
typedef pair<int,int> pii;
typedef vector<int> vi;
#define sz(n) int((n).size()) 
#define fora(cn) for(auto &x : (cn))
#define rep(i,r) for(int i=0 ;i <(r);i++)
#define rep1(i,l,r) for(int i=(l);i<(r);i++)
#define repr(i,r,l) for(int i=(r);i>=(l);i--)


// -----------------BIT -----------------
/*
 * Description: \texttt{update(i,x): a[i] += x;}\\
 * \texttt{query(i): sum in [0, i);}\\
 * \texttt{lower\_bound(sum): min pos st sum of [0, pos] >= sum, returns n if all < sum, or -1 if empty sum.
 */
struct FT {
  int n; V<ll> s;
  FT(int _n) : n(_n), s(_n) {}
  void update(int i, ll x) { 
    for (; i < n; i |= i + 1) s[i] += x; }
  ll query(int i, ll r = 0) { 
    for (; i > 0; i &= i - 1) r += s[i-1]; return r; }
  int lower_bound(ll sum) {
    if (sum <= 0) return -1; int pos = 0;
    for (int pw = 1 << (31 - __builtin_clz(n)) ; pw; pw >>= 1){
      if (pos+pw <= n && s[pos + pw-1] < sum)
        pos += pw, sum -= s[pos-1];
    }
    return pos;
  }
};


    FT f1(n), f2(n);
    // a[l...r] += v; 0 <= l <= r < n
    auto upd = [&](int l, int r, ll v) {
      f1.update(l, v), f1.update(r + 1, -v);
      f2.update(l, v*(l-1)), f2.update(r+1, -v*r);
    }; // a[l] + ... + a[r]; 0 <= l <= r < n
    auto sum = [&](int l, int r) { ++r;
      ll sub = f1.query(l) * (l-1) - f2.query(l);
      ll add = f1.query(r) * (r-1) - f2.query(r);
      return add - sub;
    };


// 2D bit with update 
template <class T=ll>
struct FT2D {
  vector<vector<T>> x;
  FT2D(int n, int m) : x(n, vector<T>(m)) { }
  void add(int k1, int k2, int a) { // x[k1][k2] += a
    for (; k1 < x.size(); k1 |= k1+1) 
      for (int k=k2; k < x[k1].size(); k |= k+1) x[k1][k] += a;
  }
  T sum(int k1, int k2 ,T s=0 ) { // return  x[0][0] to x[k1][k2] sum 
  
    for (; k1 >= 0; k1 = (k1&(k1+1))-1) 
      for (int k=k2; k >= 0; k = (k&(k+1))-1) s += x[k1][k];
    return s;
  }
  // sum of rectangle [(x1, y1), (x2, y2)] 0 base 
   T rsum(int x1, int y1, int x2, int y2) {
        return sum(x2, y2) - (x1 ? sum(x1 - 1, y2) : 0) - (y1 ? sum(x2, y1 - 1) : 0)
             + (x1 && y1 ? sum(x1 - 1, y1 - 1) : 0);
    }
};
//--------------------------------------------

// ---- data ----
struct Dsu{
   vi pa,sz;
    Dsu(int n){ 
      pa.resize(n);
      sz=vi(n,1);
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
// --------------



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
// basis
vector<int > bit(32) ;
int Sz=0 ;
void Insert(int x){
   for(int i=30;i>=0 ;i--){
       
       if((x>>i& 1) == 0 )continue ;
       
       if(bit[i] )  x^=bit[i] ;
       else {
        bit[i]=x ;
        Sz++ ;
        break ;
     }
   }
}
//------------------
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
// ------------------
V<V<ll>> PrefixSums2D(V<V<ll>> &a) {
  // ...1 base
    ll n = a.size(),m = a[0].size();
    V<V<ll>> prefixSum(n + 1, vl(m + 1, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            prefixSum[i + 1][j + 1] = prefixSum[i][j + 1] + prefixSum[i + 1][j] - prefixSum[i][j] + a[i][j];
        }
    }
    return prefixSum;
}
ll get(int rx,int ry,int lx,int ly,vvl &pre){
  // ... 1 base l to r inclusive
  lx--,--ly;
  return pre[rx][ry]+pre[lx][ly]-pre[rx][ly]-pre[lx][ry];

}

/*  pbds 
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
template<class T>using ordered_set=tree<T, null_type, less<T>, rb_tree_tag,tree_order_statistics_node_update>;
template<class T>using ordered_multiset=tree<T, null_type,less_equal<T>, rb_tree_tag,tree_order_statistics_node_update>;

*/
/// mod--- 

/**
 * Description: Finds two integers $x$ and $y$, such that $ax+by=\gcd(a,b)$. If
 * you just need gcd, use the built in \texttt{\_\_gcd} instead.
 * If $a$ and $b$ are coprime, then $x$ is the inverse of $a \pmod{b}$.
 */

ll euclid(ll a, ll b, ll &x, ll &y) {
  if (!b) return x = 1, y = 0, a;
  ll d = euclid(b, a % b, y, x);
  return y -= a/b * x, d;
}
/// Return true if there exist such (x, y) satisfy ax + by = c
/// Find (&g) = gcd(a, b)
/// Find (&x, &y) satisfy ax + by = c
// for many solution  (x0 = x0+ k*b/g ) and (y0 = y0- k* a/g)  k is  a integer 
bool find_any_solution(int a, int b, int c, int &x0, int &y0, int &g) {
    g = euclid (abs(a), abs(b), x0, y0);
    if (c % g) return false;
    x0 *= c / g;   y0 *= c / g;
    if (a < 0) x0 = -x0;
    if (b < 0) y0 = -y0;
    return true;
}
/**
 * \texttt{crt(a, m, b, n)} computes $x$ such that $x\equiv a \pmod m$, $x\equiv b \pmod n$.
 * If $|a| < m$ and $|b| < n$, $x$ will obey $0 \le x < \text{lcm}(m, n)$.
 * Assumes $mn < 2^{62}$.
 */
ll crt(ll a, ll m, ll b, ll n) {
  if (n > m) swap(a, b), swap(m, n);
  ll x, y, g = euclid(m, n, x, y);
  assert((a - b) % g == 0); // else no solution
  x = (b - a) % n * x % n / g * m + a;
  return x < 0 ? x + m*n/g : x;
}

const ll mod = 1e9+7; // change to something else
struct Mint {
  ll x;
  Mint() : x(0){}
  Mint(ll x_):x(x_ % mod){if (x < 0) x += mod;}
  Mint operator+(Mint b) { return Mint((x + b.x) % mod); }
  Mint operator-(Mint b) { return Mint((x - b.x + mod) % mod); }
  Mint operator*(Mint b) { return Mint((x * b.x) % mod); }
  Mint operator/(Mint b) { return *this * invert(b); }
  Mint invert(Mint a) {
    ll x, y, g = euclid(a.x, mod, x, y);
    assert(g == 1); return Mint((x + mod) % mod);
  }
  Mint power(ll e) {
    if (!e) return Mint(1);
    Mint r = this-> power (e / 2); r = r * r;
    return e&1 ? *this * r : r;
  }
};

vector<Mint> fact(1, 1), inv_fact(1, 1);
 
 Mint C(int n, int k) {
  if (k < 0 || k > n) return 0;
  
  while ( sz(fact) < n + 1) {
    fact.pb(fact.back() * sz(fact) );
    inv_fact.pb( Mint(1) / fact.back());
  }
  return fact[n] * inv_fact[k] * inv_fact[n - k];

}

//seive 
int preclc = 1;
vector<int> least= {0,1},primes;
void RunSieve(int n) {
  n = max(n, 1);
  least.assign(n + 1, 0);
  primes.clear();
  for (int i = 2; i <= n; i++) {
    if (least[i] == 0) {
        least[i] = i;
      primes.push_back(i);
    }
    for (int x : primes) {
      if (x > least[i] || i * x > n) {
        break;
      }
      least[i * x] = x;
    }
  }
  preclc = n;
}
vector<pii > Factorize(int  n) {
       assert (n <= preclc)  ;
       
        vector<pii > ret;
        while (n > 1) {
          if (!ret.empty() && ret.back().first == least[n]) {
            ret.back().second++;
          } else {
            ret.emplace_back(least[n], 1);
          }
          n /= least[n];
        }
    return ret;
}
/----------------------------/
// fast factor 

// works in sqrt(sqrt(n)) 
typedef unsigned long long ull;
ull modmul(ull a, ull b, ull M) {
  ll ret = a * b - M * ull(1.L / M * a * b);
  return ret + M * (ret < 0) - M * (ret >= (ll)M);
}
ull modpow(ull b, ull e, ull mod) {
  ull ans = 1;
  for (; e; b = modmul(b, b, mod), e /= 2)
    if (e & 1) ans = modmul(ans, b, mod);
  return ans;
}

bool isPrime(ull n) {
  if (n < 2 || n % 6 % 4 != 1) return (n | 1) == 3;
  ull A[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022},
      s = __builtin_ctzll(n-1), d = n >> s;
  for (ull a : A) {   // ^ count trailing zeroes
    ull p = modpow(a%n, d, n), i = s;
    while (p != 1 && p != n - 1 && a % n && i--)
      p = modmul(p, p, n);
    if (p != n-1 && i != s) return 0;
  }
  return 1;
}

ull pollard(ull n) {
  ull x = 0, y = 0, t = 30, prd = 2, i = 1, q;
  auto f = [&](ull x) { return modmul(x, x, n) + i; };
  while (t++ % 40 || __gcd(prd, n) == 1) {
    if (x == y) x = ++i, y = f(x);
    if ((q = modmul(prd, max(x,y) - min(x,y), n))) prd = q;
    x = f(x), y = f(f(y));
  }
  return gcd(prd, n);
}
vector<ull> factor(ull n) {
  if (n == 1) return {};
  if (isPrime(n)) return {n};
  ull x = pollard(n);
  auto l = factor(x), r = factor(n / x);
  l.insert(l.end(), all(r));
  return l;
}

// --------------- graph -----------

 
template<class T> struct RMQ {
  V<V<T>> jmp;
  RMQ(const V<T>& V) : jmp(1, V) {
    for (int pw = 1, k = 1; pw * 2 <= sz(V); pw *= 2, ++k) {
      jmp.emplace_back(sz(V) - pw * 2 + 1);
      rep(j,sz(jmp[k]))
        jmp[k][j] = min(jmp[k - 1][j], jmp[k - 1][j + pw]);
    }
  }
  T query(int a, int b) { // [a,b)
    assert(a < b); // or return inf if a == b
    int dep = 31 - __builtin_clz(b - a);
    return min(jmp[dep][a], jmp[dep][b - (1 << dep)]);
} };
struct LCA {
  int T = 0;
  vi time, path, ret ,tour ;
  RMQ<int> rmq;
  LCA(vector<vi>& C) : time(sz(C)), rmq((dfs(C,0,-1), ret)) {}
  void dfs(vector<vi>& C, int v, int par) {
    time[v] = T++; tour.pb(v) ;
    for (int y : C[v]) if (y != par) {
      path.pb(v), ret.pb(time[v]);
      dfs(C, y, v);
    }
    tour.pb(v) ;
  }
 
  int lca(int a, int b) {
    if (a == b) return a;
    tie(a, b) = minmax(time[a], time[b]);
    return path[rmq.query(a, b)];
  }
  //dist(a,b){return depth[a] + depth[b] - 2*depth[lca(a,b)];}
};

/**
 * Description: Calculate power of two jumps in a tree,
 * to support fast upward jumps and LCAs.
 * Assumes the root node points to itself.
 */

vector<vi> treeJump(vi& P){
  int on = 1, d = 1;
  while(on < sz(P)) on *= 2, d++;
  vector<vi> jmp(d, P);
  rep1(i,1,d) rep1(j,0,sz(P))
    jmp[i][j] = jmp[i-1][jmp[i-1][j]];
  return jmp;
}

int jmp(vector<vi>& tbl, int nod, int steps){
  rep1(i,0,sz(tbl))
    if(steps&(1<<i)) nod = tbl[i][nod];
  return nod;
}

int lca(vector<vi>& tbl, vi& depth, int a, int b) {
  if (depth[a] < depth[b]) swap(a, b);
  a = jmp(tbl, a, depth[a] - depth[b]);
  if (a == b) return a;
  for (int i = sz(tbl); i--;) {
    int c = tbl[i][a], d = tbl[i][b];
    if (c != d) a = c, b = d;
  }
  return tbl[0][a];
}



// art bri
int n , timer;
V<vi> adj; 
V<pii > bri ;
vi vis, tin,Sz, low ,art ;
void dfs(int v, int p = -1) {
    vis[v] = true;
    Sz[v]=1  ;
    tin[v] = low[v] = timer++;
    int child =0 ;
    bool pas = false;
    for (int to : adj[v]) {
        if (to == p && !pas) {
            pas = true; 
            continue;
        }
        if (vis[to]) {
            low[v] = min(low[v], tin[to]);
        } else {
            dfs(to, v);
            Sz[v]+=Sz[to] ;
            low[v] = min(low[v], low[to]);
            if (low[to] >= tin[v] && p!=-1)art[v]=1;
            if (low[to] > tin[v])bri.pb({v, to});
             ++child;
        }
    }
    if(p == -1 && child > 1)art[v ] =1 ;
}

void art_bri() {

    timer = 0; bri.clear() ;
    Sz=vi(n,0) ;
    vis=vi(n, 0);art=vi(n,0);
    tin= vi (n, -1);low=vi(n, -1);
    rep (i,n) {
      if(!vis[i])dfs(i);
    }

}

pair<vi,vi> Fun (const V<vi> &g ,const vi &order  ){
     vi topo ,vis(sz(g)),was(sz(g))  ;
     int attempt=0;
     auto dfs=[&](auto dfs,int v)-> void {
             vis[v]=1;
             was[v]=attempt ;
            for(auto u: g[v]){
               if(vis[u]==0)
                  dfs(dfs, u);
            }
        topo.pb(v);
     };
     rep(i ,sz(order) ) {
           if(vis[order[i]] ==0 ){
                attempt++ ;
              dfs(dfs, order[i] ) ;
           }
     }
    reverse(all(topo))  ;
    return make_pair(topo,was) ;
}
// return SCC 
vi SCC(const V<vi> &g ,V <vi>&rg){
    vi ord(sz(g)) ; iota(all(ord),0) ;
    auto [topo,_ ]= Fun(g, ord);
    auto [__,was]= Fun(rg,topo); // was has scc with same val 
  return was; 
}

vi topoSort(const vector<vi>& gr) {
  vi indeg(sz(gr)), q;
  for (auto& li : gr) for (int x : li) indeg[x]++;
  rep(i,sz(gr)) if (indeg[i] == 0) q.pb(i);
  rep(j,sz(q)) for (int x : gr[q[j]])
    if (--indeg[x] == 0) q.pb(x);
  return q;
}



/*
//  Floyd Warshall
  int n,m;

  vector<vector<ll>>dis(n,vector<ll>(n,inf));
  
  for(int i=0;i<n;i++)dis[i][i]=0;
    
  for(int k=0;k<n;k++){
     for(int i=0;i<n;i++){
       for(int j=0;j<n;j++){
          if(dis[i][k]!=inf && dis[k][j]!=inf ){
            dis[i][j]=min(dis[i][j],dis[i][k]+dis[k][j]);
          }
       }
     }
  }
*/

/**
 * Description: Eulerian undirected/directed path/cycle algorithm.
 * Input should be a vector of (dest, global edge index), where
 * for undirected graphs, forward/backward edges have the same index.
 * Returns a list of nodes in the Eulerian path/cycle with src at both start and end, or
 * empty list if no cycle/path exists.
 * To get edge indices back, add .second to s and ret.
 * Time: O(V + E)
 */
// adj list with edges ids 
vi eulerWalk(vector<vector<pii>>& gr, int nedges, int src=0) {
  int n = sz(gr);
  vi D(n), its(n), eu(nedges), ret, s = {src};
  D[src]++; // to allow Euler paths, not just cycles
  while (!s.empty()) {
    int x = s.back(), y, e, &it = its[x], end = sz(gr[x]);
    if (it == end){ ret.push_back(x); s.pop_back(); continue; }
    tie(y, e) = gr[x][it++];
    if (!eu[e]) {
      D[x]--, D[y]++;
      eu[e] = 1; s.push_back(y);
    }}
  for (int x : D) if (x < 0 || sz(ret) != nedges+1) return {};
  return {ret.rbegin(), ret.rend()};
}

/// ----geo--- 

template <class T> int sgn(T x) { return (x > 0) - (x < 0); }
const double eps = 1e-9;
const double PI = acos((double)-1.0);
// int sgn(double x) { return (x > eps) - (x < -eps); }

template<class T> struct Point { 
typedef Point P ;T x, y;
explicit Point(T _x=0, T _y=0) : x(_x),y(_y){}
bool operator<(P p) const { return tie(x,y) < tie(p.x,p.y); }
bool operator==(P p) const { return tie(x,y)==tie(p.x,p.y); }
P operator+(P p) const{return P(x+p.x,y+p.y);}
P operator-(P p) const{return P(x-p.x,y-p.y);}
P operator*(T d) const { return P(x*d, y*d); }
P operator/(T d) const { return P(x/d, y/d); }
T dot(P p) const { return x*p.x + y*p.y; }
T cross(P p) const { return x*p.y - y*p.x; }
T cross(P a, P b) const { return (a-*this).cross(b-*this); }
T dist2() const { return x*x + y*y; }
double dist() const { return sqrt((double)dist2()); }
// angle to x-axis in interval [-pi, pi]
double angle() const { return atan2(y, x); }
// makes dist() = 1
P unit() const { return *this/dist(); } 
// rotate by +90 degree
P perp() const { return P(-y, x); }
P normal() const { return perp().unit(); }
//rotate 'a' radians ccw around (0,0)
P rotate(double a) const { return P(x*cos(a)-y*sin(a),x*sin(a)+y*cos(a)); }
friend ostream& operator<<(ostream& os, P p) { 
  return os<<"("<< p.x << "," << p.y << ")";}
friend istream &operator>>(istream &in,P &p){return in >>p.x>>p.y; }
};
// a--b seg [cl,0] [-1, cw ] [1 ,ccw]
template<class P> 
int orientation(P a, P b, P c) { return sgn(a.cross(b,c)); }

double rad_to_deg(double r) { return (r * 180.0 / PI); }
double deg_to_rad(double d) { return (d * PI / 180.0); }

template<class P>  double get_angle(P a, P b) {
    double costheta = a.dot(b) / a.dist() / b.dist();
    return acos(max((double)-1.0, min((double)1.0, costheta)));
}
template<class P> 
bool is_point_in_angle(P b, P  a, P c, P  p) { // does point p lie in angle <bac
    assert(orientation(a, b, c) != 0);
    if (orientation(a, c, b) < 0) swap(b, c);
    return orientation(a, c, p) >= 0 && orientation(a, b, p) <= 0;
}
// checks is point p in s-- e line segment 
template<class P> 
bool onSegment(P s, P e, P p) {return p.cross(s, e) == 0 && (s - p).dot(e - p) <= 0;}
template<class P> 
vector<P> segInter(P a, P b, P c, P d) {
     auto oa = c.cross(d, a), ob = c.cross(d, b),
          oc = a.cross(b, c), od = a.cross(b, d);
     // Checks if intersection is single non-endpoint point.
     if (sgn(oa) * sgn(ob) < 0 && sgn(oc) * sgn(od) < 0)
          return {(a * ob - b * oa) / (ob - oa)};
     set<P> s;
     if (onSegment(c, d, a)) s.insert(a);
     if (onSegment(c, d, b)) s.insert(b);
     if (onSegment(a, b, c)) s.insert(c);
     if (onSegment(a, b, d)) s.insert(d);
     return {all(s)};
}
//Returns the shortest distance between point p and the line segment from point s to e.
// typedef Point<double> P;
// double segDist(P& s, P& e, P& p) {
//   if (s==e) return (p-s).dist();
//   auto d = (e-s).dist2(), t = min(d,max(.0,(p-s).dot(e-s)));
//   return ((p-s)*d-(e-s)*t).dist()/d;
// }


// [-1 for boundary ][1 for inside ] else outside
template<class P> 
int P_i_polygon(vector<P>&p ,P &x){
    P y = P(x.x + 1, ll(1e9)+7ll); 
    int cnt =0 ,n=sz(p) ;
    rep(i,n){
            cnt+=sz(segInter(p[i] ,p[(i+1)%n] ,x,y) )>0 ;
            if(onSegment(p[i], p[ (i+1)%n], x) )return -1 ;
     }
      if(cnt&1) return 1;
      else return 0 ;
}


template< class P>
vector<P> convex_hull(vector<P> &p) {
  if (p.size() <= 1) return p;
  vector<P> v = p;
    sort(v.begin(), v.end());
    vector<P> up, dn;
    for (auto& p : v) {
        while (up.size() > 1 && orientation(up[up.size() - 2], up.back(), p) >= 0) {
            up.pop_back();
        }
        while (dn.size() > 1 && orientation(dn[dn.size() - 2], dn.back(), p) <= 0) {
            dn.pop_back();
        }
        up.pb(p); dn.pb(p);
    }
    v = dn;
    if (v.size() > 1) v.pop_back();
    reverse(up.begin(), up.end());
    up.pop_back();
    for (auto& p : up)  v.pb(p);
    
    if (v.size() == 2 && v[0] == v[1]) v.pop_back();
    return v;
}



// -string -- 
//  Description: Static hashing for 0-indexed string. Intervals are $[l, r]$.
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



// Description: Hashing with point updates on string (0-indexed). 
// \texttt{upd(i, x): s[i] += x}. Intervals are $[l,r]$.
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
 


/* * Description: pi[x] computes the length of the longest prefix of s that ends at x,
 * other than s[0...x] itself (abacaba -> 0010123).
 * Can be used to find all occurrences of a string. */
vi pi(const string& s) {
  vi p(sz(s));
  rep1(i,1,sz(s)) {
    int g = p[i-1];
    while (g && s[i] != s[g]) g = p[g-1];
    p[i] = g + (s[i] == s[g]);
  }
  return p;
}

vi match(const string& s, const string& pat) {
  vi p = pi(pat + '\0' + s), res;
  rep1(i,sz(p)-sz(s),sz(p))
    if (p[i] == sz(pat)) res.push_back(i - 2 * sz(pat));
  return res;
}

/*
 * Description: z[i] computes the length of the longest common prefix of s[i:] and s,
 * except z[0] = 0. (abacaba -> 0010301)
 */
vi Z(const string& S) {
  vi z(sz(S));
  int l = -1, r = -1;
  rep1(i,1,sz(S)) {
    z[i] = i >= r ? 0 : min(r - i, z[i - l]);
    while (i + z[i] < sz(S) && S[i + z[i]] == S[z[i]])
      z[i]++;
    if (i + z[i] > r)
      l = i, r = i + z[i];
  }
  return z;
}

/**
 * Description: Finds the lexicographically smallest rotation of a string.
 * Time: O(N)
 * Usage:
 *  rotate(v.begin(), v.begin()+minRotation(v), v.end());
 */
int minRotation(string s) {
  int a=0, N=sz(s); s += s;
  rep1(b,0,N) rep1(k,0,N) {
    if (a+k == b || s[a+k] < s[b+k]) {b += max(0, k-1); break;}
    if (s[a+k] > s[b+k]) { a = b; break; }
  }
  return a;
}

/**
 * Description: For each position in a string, computes p[0][i] = half length of
 *  longest even palindrome around pos i, p[1][i] = longest odd (half rounded down).
 * Time: O(N)
 */
array<vi, 2> manacher(const string& s) {
  int n = sz(s);
  array<vi,2> p = {vi(n+1), vi(n)};
  rep1(z,0,2) for (int i=0,l=0,r=0; i < n; i++) {
    int t = r-i+!z;
    if (i<r) p[z][i] = min(t, p[z][l+t]);
    int L = i-p[z][i], R = i+p[z][i]-!z;
    while (L>=1 && R+1<n && s[L-1] == s[R+1])
      p[z][i]++, L--, R++;
    if (R>r) l=L, r=R;
  }
  return p;
}

// segtree 
template <class T, T (*op)(T, T), T (*e)()> 
struct segtree {
  int n,size;
  vector<T> d;
  segtree(const vector<T>& v){
      n=int(v.size()),size=1; 
      while(size<n)size<<=1;
      
        d=vector<T>(2 * size, e());
        for(int i=0;i<n;i++)d[size + i] = v[i];
        for(int i=size-1;i>=1;i--) update(i);
    }
    void update(int k) { 
      d[k] = op(d[2 * k], d[2 * k + 1]); 
    }
    void set(int p, T x) {
        assert(0 <= p && p < n);
        p += size;  d[p] = x;
        while(p > 1)update(p >>=1);
    }
    T prod(int l, int r) const { // [l,r)
      assert(0 <= l && l <= r && r <= n);
        T sml = e(), smr = e();
        for(l += size ,r += size ; l<r ;  l>>=1,r>>=1 ){
            if (l & 1) sml = op(sml, d[l++]);
            if (r & 1) smr = op(d[--r], smr); 
        }
      return op(sml, smr);
    }
};

ll op(ll x ,ll y){ return x+y;}
ll e(){return 0ll ;}


template <class S,
          S (*op)(S, S),
          S (*e)(),
          class F,
          S (*Map)(F, S),
          F (*Com)(F, F),
          F (*id)()>
struct lazy_segtree {
 
    lazy_segtree(const std::vector<S>& v) : _n(int(v.size())) {  
        size = 1 ;while(size<_n)size*=2 ;
        log = 31 - __builtin_clz(size); 
        d = std::vector<S>(2 * size, e());
        lz = std::vector<F>(size, id());
        for (int i = 0; i < _n; i++) d[size + i] = v[i];
        for (int i = size - 1; i >= 1; i--) update(i);
        
    }
 
    void set(int p, S x) {
        assert(0 <= p && p < _n);
        p += size;
        for (int i = log; i >= 1; i--) push(p >> i);
        d[p] = x;
        for (int i = 1; i <= log; i++) update(p >> i);
    }
 
    S get(int p) {
        assert(0 <= p && p < _n);
        p += size;
        for (int i = log; i >= 1; i--) push(p >> i);
        return d[p];
    }
 
    S prod(int l, int r) { // [l,r)
        assert(0 <= l && l <= r && r <= _n);
        if (l == r) return e();
 
        l += size;
        r += size;
 
        for (int i = log; i >= 1; i--) {
            if (((l >> i) << i) != l) push(l >> i);
            if (((r >> i) << i) != r) push((r - 1) >> i);
        }
 
        S sml = e(), smr = e();
        while (l < r) {
            if (l & 1) sml = op(sml, d[l++]);
            if (r & 1) smr = op(d[--r], smr);
            l >>= 1;
            r >>= 1;
        }
 
        return op(sml, smr);
    }
 
    S all_prod() { return d[1]; }
 
    void apply(int p, F f) {
        assert(0 <= p && p < _n);
        p += size;
        for (int i = log; i >= 1; i--) push(p >> i);
        d[p] = Map(f, d[p]);
        for (int i = 1; i <= log; i++) update(p >> i);
    }
    void apply(int l, int r, F f) { // [l,r)
        assert(0 <= l && l <= r && r <= _n);
        if (l == r) return;
 
        l += size;
        r += size;
 
        for (int i = log; i >= 1; i--) {
            if (((l >> i) << i) != l) push(l >> i);
            if (((r >> i) << i) != r) push((r - 1) >> i);
        }
 
        {
            int l2 = l, r2 = r;
            while (l < r) {
                if (l & 1) all_apply(l++, f);
                if (r & 1) all_apply(--r, f);
                l >>= 1;
                r >>= 1;
            }
            l = l2;
            r = r2;
        }
 
        for (int i = 1; i <= log; i++) {
            if (((l >> i) << i) != l) update(l >> i);
            if (((r >> i) << i) != r) update((r - 1) >> i);
        }
    }
 
    template <bool (*g)(S)> int max_right(int l) {
        return max_right(l, [](S x) { return g(x); });
    }
    template <class G> int max_right(int l, G g) {
        assert(0 <= l && l <= _n);
        assert(g(e()));
        if (l == _n) return _n;
        l += size;
        for (int i = log; i >= 1; i--) push(l >> i);
        S sm = e();
        do {
            while (l % 2 == 0) l >>= 1;
            if (!g(op(sm, d[l]))) {
                while (l < size) {
                    push(l);
                    l = (2 * l);
                    if (g(op(sm, d[l]))) {
                        sm = op(sm, d[l]);
                        l++;
                    }
                }
                return l - size;
            }
            sm = op(sm, d[l]);
            l++;
        } while ((l & -l) != l);
        return _n;
    }
 
    template <bool (*g)(S)> int min_left(int r) {
        return min_left(r, [](S x) { return g(x); });
    }
    template <class G> int min_left(int r, G g) {
        assert(0 <= r && r <= _n);
        assert(g(e()));
        if (r == 0) return 0;
        r += size;
        for (int i = log; i >= 1; i--) push((r - 1) >> i);
        S sm = e();
        do {
            r--;
            while (r > 1 && (r % 2)) r >>= 1;
            if (!g(op(d[r], sm))) {
                while (r < size) {
                    push(r);
                    r = (2 * r + 1);
                    if (g(op(d[r], sm))) {
                        sm = op(d[r], sm);
                        r--;
                    }
                }
                return r + 1 - size;
            }
            sm = op(d[r], sm);
        } while ((r & -r) != r);
        return 0;
    }
 
  private:
    int _n, size, log;
    std::vector<S> d;
    std::vector<F> lz;
 
    void update(int k) { d[k] = op(d[2 * k], d[2 * k + 1]); }
    void all_apply(int k, F f) {
        d[k] = Map(f, d[k]);
        if (k < size) lz[k] = Com(f, lz[k]);
    }
    void push(int k) {
        all_apply(2 * k, lz[k]);
        all_apply(2 * k + 1, lz[k]);
        lz[k] = id();
    }
};
 
// [[ Assignment and Sum ]]   
using S = pll  ; using F = ll ;
S op( S a,S b) {return {a.ff+b.ff,a.ss+b.ss} ;}
S e() {return  {0,0}; }
S Map(F f, S x) { if(f==-1)return  x;
  else return {f*x.ss ,x.ss};}
F Com(F f, F  g) { return  f==-1?g:f ;}
F id() { return -1 ; }

// [[ Addition and Sum ]] 
using S = pll; using F = ll ;
S op( S a,S b) { return  {a.ff+b.ff,a.ss+b.ss};}
S e() { return  {0,0}; }
S Map(F f, S x) {return {x.ff + (x.ss*f ),x.ss};}
F Com(F f, F g) { return  f+g ; }
F id() { return 0; }

// psegtree 
struct node{
  int64_t sum ;
  node *lc; node *rc;
  node(){ rc=lc=NULL;sum =0 ;
  }
};
node* build(int L,int R){
  node * cur  = new node() ;      
   if(L == R){
      cur->sum=a[L ];
     return  cur;
   }
   int m=(L+R)/2;
   cur->lc = build(L,m);
   cur->rc= build(m+1,R);
   cur->sum = cur->lc-> sum + cur->rc->sum ;  
   return cur ;  
}
ll quary (node *cur, int a,int b, int L,int R){
     
    if(a>R  || b< L )return 0ll ;
    if(L>= a && R<=b )return cur-> sum ;
   int m=(L+R)/2;
   ll p1 =quary(cur->lc ,a,b ,L,m);
   ll q1 =quary( cur->rc,a,b ,m+1,R);
  return p1+ q1 ;
}
node* update(node *cur ,int i,int val ,int L,int R){
   if(L==R ){
     node * ne = new node() ;
     ne ->sum= val ;
     return  ne ;
   }
   int m=(L+R)/2;
   if( i<= m ){
    node * ne  = new node() ;
    ne-> lc =update(cur->lc, i,val ,L,m) ;
    ne-> rc = cur-> rc ;
    ne-> sum = ne->lc-> sum + ne->rc->sum ;  
    return ne ;
  }else{
    node * ne  = new node() ;
    ne-> rc =update(cur->rc, i,val ,m+1,R ) ;
  
    ne-> lc = cur-> lc ;
    ne-> sum = ne->lc-> sum + ne->rc->sum ;  
     return ne ; 
   }
  
}

// ----HLd-----
// build  ord[ Id[i] ] =  a[i] for  tree 
template <bool VALS_EDGES = false >
struct HLD {
  int N, tim = 0;
  vector<vi> adj;
  vi par, siz, rt, pos;

  HLD(vector<vi> adj_)
    : N(sz(adj_)), adj(adj_), par(N, -1), siz(N, 1),
      rt(N),pos(N) { dfsSz(0); dfsHld(0); }
  void dfsSz(int v) {
    for (int& u : adj[v]) {
      adj[u].erase(find(all(adj[u]), v));
      par[u] = v;
      dfsSz(u);
      siz[v] += siz[u];
      if (siz[u] > siz[adj[v][0]]) swap(u, adj[v][0]);
    }
  }
  void dfsHld(int v) {
    pos[v] = tim++;
    for (int u : adj[v]) {
      rt[u] = (u == adj[v][0] ? rt[v] : u);
      dfsHld(u);
    }
  }
  void process(int u, int v, vii &op) {
    for (;; v = par[rt[v]]) {
      if (pos[u] > pos[v]) swap(u, v);
      if (rt[u] == rt[v]) break;
      op.pb({pos[rt[v]], pos[v] + 1});
    }
     op.pb({pos[u] + VALS_EDGES, pos[v] + 1});
  }
  vii Path(int u, int v) { // return vec of [l,r)  
   vii op ; process(u, v,op) ;
    return op ;
  }
  int Id (int v){ //  return i such ord[i ] = a[i]  for tree 
     return pos[v] ;
  }
  pii querySubtree(int v) { // return [l,r)    
    return  {pos[v] + VALS_EDGES, pos[v] + siz[v]} ;
  }
};

// mo 

void MO_qary(){
   
/*
  mo tree quary  make 2n tour vec 
  for path u to v if  u == lca(u,v) then quary = in[u] to in [v ] ;
  else quary  out[u] to in[v] then add lca node  
  all are 0 base 
*/

const int Block = 335  ;
const int N =  100009;

struct  MO{
  int l,r , id ;
  bool operator<(const MO &e )const {
   return  make_pair(l/Block , r)< make_pair(e.l/Block ,e.r) ;
  }
} quary[N];
/*
void add (int id){   }
void rem(int id){ }
int clc(){ }
*/
vi tour , in(n),out(n) ;
 
    auto dfs=[&](auto dfs,int v,int pa)-> void {
        in[v]=sz(tour);
        tour.pb(v) ;
        for(auto u: G[v]){
           if(u!= pa ){
             dfs(dfs,u,v) ;
           }
        }
       out[v]=sz(tour);
       tour.pb(v);
    }; dfs(dfs,0,-1);
    
  

 rep(i,q ){ // number of quary 
  int u ,v ;
   cin>>u>>v ;
   --u,--v;
    if(in[u]> in[v] )swap(u,v) ;
      int x= A.lca(u,v) ;
    if( x == u){
       quary[i]={in[u],in[v],i};
    }else{
       quary[i]={out[u],in[v],i};
    }
 }
 sort(quary,quary+q);
  
int cl  = 0 ,cr=  -1 ;
vi ans (q) ;
for (int i=0;i<q;i++) {
  auto &[l,r ,k, lca, id ]  = quary[i];     
    while (cl > l)add( --cl);
    while (cr < r)add(++cr );
    while (cl < l)rem(cl++ );
    while (cr > r)rem (cr-- );
    ans[id]= clc () ;
}


}
 

void MO_qary(){
  // mo  quary 0 base
const int Block = 335 ,N =  100009;

struct  MO{
   int l,r , id ;
   bool operator<(const MO &e )const {
     return  make_pair(l/Block , r)< make_pair(e.l/Block ,e.r) ;
    }
} quary[N];
/*
void add (int id){   }
void rem(int id){ }
int clc(){ }
*/
rep(i,q ){ // number of quary 
 
   int l,r ;
   cin>>l>>r;
   --l,--r; 
   quary[i]={l,r,i};
 }
 sort(quary,quary+q);
 int cl  = 0 ,cr=  -1 ;
vi ans (q) ;
for (int i=0;i<q;i++) {
  auto &[l,r ,k, lca, id ]  = quary[i];     
    while (cl > l)add( --cl);
    while (cr < r)add(++cr );
    while (cl < l)rem(cl++ );
    while (cr > r)rem (cr-- );
    ans[id]= clc () ;
}


} 
// segtree sort 

template<class S> struct segtree {
  int n; V<S> t;
  void init(int _) { n = _; t.assign(n+n-1, S()); }
  void init(const V<S>& v) { 
    n = sz(v); t.assign(n + n - 1, S()); 
    build(0,0,n-1,v); 
  } template <typename... T>
  void upd(int l, int r, const T&... v) {
    assert(0 <= l && l <= r && r < n);
    upd(0, 0, n-1, l, r, v...);
  }
  S get(int l, int r) {
    assert(0 <= l && l <= r && r < n);
    return get(0, 0, n-1, l, r);
  }
private:
  inline void push(int u, int b, int e) {
    if (t[u].lazy == 0) return;
    int mid = (b+e)>>1, rc = u+((mid-b+1)<<1);
    t[u+1].upd(b, mid, t[u].lazy);
    t[rc].upd(mid+1, e, t[u].lazy);
    t[u].lazy = 0;
  }
  void build(int u,int b,int e,const V<S>&v) {
    if (b == e) return void(t[u] = v[b]);
    int mid = (b+e)>>1, rc = u+((mid-b+1)<<1);
    build(u+1, b,mid,v); build(rc, mid+1,e,v);
    t[u] = t[u+1] + t[rc];
  } template<typename... T>
  void upd(int u, int b, int e, int l, int r, const T&... v) {
    if (l <= b && e <= r) return t[u].upd(b, e, v...);
    push(u, b, e);
    int mid = (b+e)>>1, rc = u+((mid-b+1)<<1);
    if (l<=mid) upd(u+1, b, mid, l, r, v...);
    if (mid<r) upd(rc, mid+1, e, l, r, v...);
    t[u] = t[u+1] + t[rc];
  }
  S get(int u, int b, int e, int l, int r) {
    if (l <= b && e <= r) return t[u];
    push(u, b, e); 
    S res; int mid = (b+e)>>1, rc = u+((mid-b+1)<<1);
    if (r<=mid) res = get(u+1, b, mid, l, r);
    else if (mid<l) res = get(rc,mid+1,e,l,r);
    else res = get(u+1, b, mid, l, r) + get(rc, mid+1, e, l, r);
    t[u] = t[u+1] + t[rc]; return res;
  }
}; // Hash upto here = 773c09
/* (1) Declaration:
Create a node class. Now, segtree<node> T;
T.init(10) creates everything as node()
Consider using V<node> leaves to build
(2) upd(l, r, ...v): update range [l, r]
order in ...v must be same as node.upd() fn */
struct node {
  ll sum = 0, lazy = 0;
  node () {} // write full constructor
  node operator+(const node &obj) { 
    return {sum + obj.sum, 0};    }
  void upd(int b, int e, ll x) {
    sum += (e - b + 1) * x, lazy += x;
} }; 
