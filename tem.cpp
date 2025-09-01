/*

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

*/
#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <chrono>
#include <climits>
#include <cstdint>
#include <cmath>
#include <complex>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <stack>
#include <queue>
#include <random>
#include <set>
#include <vector>
using namespace std;
 


template<typename A, typename B> ostream& operator<<(ostream &os, const pair<A, B> &p) { return os << '(' << p.first << ", " << p.second << ')'; }
template<typename... Args> ostream& operator<<(ostream& os, const tuple<Args...>& t) { os << '('; apply([&os](const Args&... args) { size_t n = 0; ((os << args << (++n != sizeof...(Args) ? ", " : "")), ...); }, t); return os << ')'; }
template<typename T_container, typename T = typename enable_if<!is_same<T_container, string>::value, typename T_container::value_type>::type> ostream& operator<<(ostream &os, const T_container &v) { os << '{'; string sep; for (const T &x : v) os << sep << x, sep = ", "; return os << '}'; }
void dbg_out() { cerr << endl; }
template<typename Head, typename... Tail> void dbg_out(Head H, Tail... T) { cerr << ' ' << H; dbg_out(T...); }
#ifndef ONLINE_JUDGE 
#define dbg(...) cerr << '[' << __LINE__ << "] (" << #__VA_ARGS__ << "):", dbg_out(__VA_ARGS__)
#else
#define dbg(...)
#endif
 
#define Tct   template<typename T>
#define Tctu  template<typename T,typename U>
#define clr(cnt, x) memset((cnt), (x), sizeof(cnt))
typedef long long ll;
Tct using V=vector<T>;
#define all(x) (x).begin(), (x).end()
#define lb lower_bound
#define ub upper_bound
#define bg(cn) cn.begin()
#define Pop(cn) cn.pop_back()
#define ff first
#define bk(x) x.back() 
#define fk(x) x.fornt()
#define eb emplace_back
#define ss second 
#define pb push_back
 
Tct void remDup(vector<T> &v){sort(all(v)); v.erase(unique(all(v)),end(v));}
Tctu void safeErase(T &t, const U &u) {auto it = t.find(u);assert(it != end(t)); t.erase(it);}
Tct int lwb(const V<T>&a,const T&b){return int(lb(all(a), b)-bg(a));}// (v[r]>=x && v[l]<x) return r;
Tct int upb(const V<T>&a,const T&b){return int(ub(all(a),b)-bg(a));}//( v[r]>x && v[l]<=x) return r;
Tctu int CC(const T&a,const U&b){return int(count(all(a),b));}
mt19937 rng(chrono::high_resolution_clock::now().time_since_epoch().count());
inline ll getrandom(ll a,ll b) { return uniform_int_distribution<ll>(a,b)(rng); }
 
void YN(bool x){if(x)cout<<"YES\n";else cout<<"NO\n";}
ll SUM_(ll n){return (n*(n+1))/2;}
ll SUM_(ll a,ll b) { return (b- a+1) * (a + b) / 2;}
ll A_SUM(ll n,ll a,ll b){return (n*(a+b))/2;}
ll G_SUM(ll k,ll a,ll b){return ((b*k)-a)/(k-1);}
ll cdiv(ll a, ll b) {return a / b + ((a ^ b) > 0 && a % b);} 
Tctu bool cmax(T &a, U b){return (a<b?a=b,1:0 );}
Tctu bool cmin(T &a, U b){return (a>b?a=b,1:0 );}
const ll inf = 2e18;
#define PQ priority_queue
Tct using  PQ1 =PQ<T,V<T>,greater<T>>;
 
#define to_str(x) (to_string(x)) 
#define to_ll(x) (stoll(x))
#define to_char(x) ((x)+'0')
#define to_int(x) ((x)-'0')
 
const int dx[] = {1, -1, 0, 0, 1, 1, -1, -1};  
const int dy[] = {0, 0, 1, -1, 1, -1, 1, -1};  
 
constexpr int Pct(ll x){return __builtin_popcountll(x);}  
constexpr int Msb(ll x){return int(log2(x));}  
constexpr __int128_t P2(ll x) { return __int128_t(1) << x; }
constexpr ll Msk2(ll x) { return P2(x) - 1;}
constexpr ll Lsb(ll x) { return  x&~(x-1);}
constexpr  bool Is_set(ll x,ll j) { return  (x&(1ll<<j) )==(1ll<<j);}
 
typedef pair<int,int> pii;
typedef pair<ll, ll> pll;
typedef V<pll> vll;
typedef V<pii >vii;
typedef V<V<int> > vvi;
typedef vector<ll>vl;  
typedef V<vl>vvl;
#define sz(n) int((n).size()) 
#define fora(cn) for(auto &x : (cn))
using   u_128 = unsigned __int128;
using   i_128=__int128_t;
#define PI 2*acos(0) 
#define el '\n'
typedef unsigned long long ull;
typedef vector<int> vi;
#define se(X) setprecision(X)
#define rall(x) (x).rbegin(), (x).rend()
#define rep(i,r) for(int (i)=(0);(i)<(r);(i)++)
#define rep1(i,l,r) for(int (i)=(l);(i)<(r);(i)++)
#define repr(i,r,l) for(int (i)=(r);(i)>=(l);(i)--)
const int MX=100005;

void solve(){
  
}
int main() {
 ios::sync_with_stdio(false);
 cout.tie(nullptr) ;
 cin.tie(nullptr);
 
 int tt=1;
   cin>>tt;
  while(tt--){
    solve();
  }
 return 0;
}   
 
 

