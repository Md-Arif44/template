#include<bits/stdc++.h>
using namespace std;
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
#ifndef DEBUG_TEMPLATE_CPP
#define DEBUG_TEMPLATE_CPP
namespace __DEBUG_UTIL__
{
    using namespace std;
    /* Primitive Datatypes Print */
    void print(const char *x) { cerr << x; }
    void print(bool x) { cerr << (x ? "T" : "F"); }
    void print(char x) { cerr << '\'' << x << '\''; }
    void print(signed short int x) { cerr << x; }
    void print(unsigned short int x) { cerr << x; }
    void print(signed int x) { cerr << x; }
    void print(unsigned int x) { cerr << x; }
    void print(signed long int x) { cerr << x; }
    void print(unsigned long int x) { cerr << x; }
    void print(signed long long int x) { cerr << x; }
    void print(unsigned long long int x) { cerr << x; }
    void print(float x) { cerr << x; }
    void print(double x) { cerr << x; }
    void print(long double x) { cerr << x; }
    void print(string x) { cerr << '\"' << x << '\"'; }
    template <size_t N>
    void print(bitset<N> x) { cerr << x; }
    void print(vector<bool> v)
    { /* Overloaded this because stl optimizes vector<bool> by using
          _Bit_reference instead of bool to conserve space. */
        int f = 0;
        cerr << '{';
        for (auto &&i : v)
            cerr << (f++ ? "," : "") << (i ? "T" : "F");
        cerr << "}";
    }
    /* Templates Declarations to support nested datatypes */
    template <typename T>
    void print(T &&x);
    template <typename T>
    void print(vector<vector<T>> mat);
    template <typename T, size_t N, size_t M>
    void print(T (&mat)[N][M]);
    template <typename F, typename S>
    void print(pair<F, S> x);
    template <typename T, size_t N>
    struct Tuple;
    template <typename T>
    struct Tuple<T, 1>;
    template <typename... Args>
    void print(tuple<Args...> t);
    template <typename... T>
    void print(priority_queue<T...> pq);
    template <typename T>
    void print(stack<T> st);
    template <typename T>
    void print(queue<T> q);
    /* Template Datatypes Definitions */
    template <typename T>
    void print(T &&x)
    {
        /*  This works for every container that supports range-based loop
            i.e. vector, set, map, oset, omap, dequeue */
        int f = 0;
        cerr << '{';
        for (auto &&i : x)
            cerr << (f++ ? "," : ""), print(i);
        cerr << "}";
    }
    template <typename T>
    void print(vector<vector<T>> mat)
    {
        int f = 0;
        cerr << "\n~~~~~\n";
        for (auto &&i : mat)
        {
            cerr << setw(2) << left << f++, print(i), cerr << "\n";
        }
        cerr << "~~~~~\n";
    }
    template <typename T, size_t N, size_t M>
    void print(T (&mat)[N][M])
    {
        int f = 0;
        cerr << "\n~~~~~\n";
        for (auto &&i : mat)
        {
            cerr << setw(2) << left << f++, print(i), cerr << "\n";
        }
        cerr << "~~~~~\n";
    }
    template <typename F, typename S>
    void print(pair<F, S> x)
    {
        cerr << '(';
        print(x.first);
        cerr << ',';
        print(x.second);
        cerr << ')';
    }
    template <typename T, size_t N>
    struct Tuple
    {
        static void printTuple(T t)
        {
            Tuple<T, N - 1>::printTuple(t);
            cerr << ",", print(get<N - 1>(t));
        }
    };
    template <typename T>
    struct Tuple<T, 1>
    {
        static void printTuple(T t) { print(get<0>(t)); }
    };
    template <typename... Args>
    void print(tuple<Args...> t)
    {
        cerr << "(";
        Tuple<decltype(t), sizeof...(Args)>::printTuple(t);
        cerr << ")";
    }
    template <typename... T>
    void print(priority_queue<T...> pq)
    {
        int f = 0;
        cerr << '{';
        while (!pq.empty())
            cerr << (f++ ? "," : ""), print(pq.top()), pq.pop();
        cerr << "}";
    }
    template <typename T>
    void print(stack<T> st)
    {
        int f = 0;
        cerr << '{';
        while (!st.empty())
            cerr << (f++ ? "," : ""), print(st.top()), st.pop();
        cerr << "}";
    }
    template <typename T>
    void print(queue<T> q)
    {
        int f = 0;
        cerr << '{';
        while (!q.empty())
            cerr << (f++ ? "," : ""), print(q.front()), q.pop();
        cerr << "}";
    }
    /* Printer functions */
    void printer(const char *) {} /* Base Recursive */
    template <typename T, typename... V>
    void printer(const char *names, T &&head, V &&...tail)
    {
        /* Using && to capture both lvalues and rvalues */
        int i = 0;
        for (size_t bracket = 0; names[i] != '\0' and (names[i] != ',' or bracket != 0); i++)
            if (names[i] == '(' or names[i] == '<' or names[i] == '{')
                bracket++;
            else if (names[i] == ')' or names[i] == '>' or names[i] == '}')
                bracket--;
        cerr.write(names, i) << " = ";
        print(head);
        if (sizeof...(tail))
            cerr << " ||", printer(names + i + 1, tail...);
        else
            cerr << "]\n";
    }
    /* PrinterArr */
    void printerArr(const char *) {} /* Base Recursive */
    template <typename T, typename... V>
    void printerArr(const char *names, T arr[], size_t N, V... tail)
    {
        size_t ind = 0;
        for (; names[ind] and names[ind] != ','; ind++)
            cerr << names[ind];
        for (ind++; names[ind] and names[ind] != ','; ind++)
            ;
        cerr << " = {";
        for (size_t i = 0; i < N; i++)
            cerr << (i ? "," : ""), print(arr[i]);
        cerr << "}";
        if (sizeof...(tail))
            cerr << " ||", printerArr(names + ind + 1, tail...);
        else
            cerr << "]\n";
    }
}
#ifndef ONLINE_JUDGE
#define debug(...) std::cerr << __LINE__ << ": [", __DEBUG_UTIL__::printer(#__VA_ARGS__, __VA_ARGS__)
#define debugArr(...) std::cerr << __LINE__ << ": [", __DEBUG_UTIL__::printerArr(#__VA_ARGS__, __VA_ARGS__)
#else
#define debug(...)
#define debugArr(...)
#endif
#endif
using namespace __gnu_pbds;
template<class T>using ordered_set=tree<T, null_type, less<T>, rb_tree_tag,tree_order_statistics_node_update>;
template<class T>using ordered_multiset=tree<T, null_type,less_equal<T>, rb_tree_tag,tree_order_statistics_node_update>;
typedef long long ll;
// template <typename T>
// T inverse(T a, T m) {
//   T u = 0, v = 1;
//   while (a != 0) {
//     T t = m / a;
//     m -= t * a; swap(a, m);
//     u -= t * v; swap(u, v);
//   }
//   assert(m == 1);
//   return u;
// }

// template <typename T>
// class Modular {
//  public:
//   using Type = typename decay<decltype(T::value)>::type;

//   constexpr Modular() : value() {}
//   template <typename U>
//   Modular(const U& x) {
//     value = normalize(x);
//   }

//   template <typename U>
//   static Type normalize(const U& x) {
//     Type v;
//     if (-mod() <= x && x < mod()) v = static_cast<Type>(x);
//     else v = static_cast<Type>(x % mod());
//     if (v < 0) v += mod();
//     return v;
//   }

//   const Type& operator()() const { return value; }
//   template <typename U>
//   explicit operator U() const { return static_cast<U>(value); }
//   constexpr static Type mod() { return T::value; }

//   Modular& operator+=(const Modular& other) { if ((value += other.value) >= mod()) value -= mod(); return *this; }
//   Modular& operator-=(const Modular& other) { if ((value -= other.value) < 0) value += mod(); return *this; }
//   template <typename U> Modular& operator+=(const U& other) { return *this += Modular(other); }
//   template <typename U> Modular& operator-=(const U& other) { return *this -= Modular(other); }
//   Modular& operator++() { return *this += 1; }
//   Modular& operator--() { return *this -= 1; }
//   Modular operator++(int) { Modular result(*this); *this += 1; return result; }
//   Modular operator--(int) { Modular result(*this); *this -= 1; return result; }
//   Modular operator-() const { return Modular(-value); }

//   template <typename U = T>
//   typename enable_if<is_same<typename Modular<U>::Type, int>::value, Modular>::type& operator*=(const Modular& rhs) {
//     value = normalize(static_cast<int64_t>(value) * static_cast<int64_t>(rhs.value));
//     return *this;
//   }
//   template <typename U = T>
//   typename enable_if<is_same<typename Modular<U>::Type, long long>::value, Modular>::type& operator*=(const Modular& rhs) {
//     long long q = static_cast<long long>(static_cast<long double>(value) * rhs.value / mod());
//     value = normalize(value * rhs.value - q * mod());
//     return *this;
//   }
//   template <typename U = T>
//   typename enable_if<!is_integral<typename Modular<U>::Type>::value, Modular>::type& operator*=(const Modular& rhs) {
//     value = normalize(value * rhs.value);
//     return *this;
//   }

//   Modular& operator/=(const Modular& other) { return *this *= Modular(inverse(other.value, mod())); }

//   friend const Type& abs(const Modular& x) { return x.value; }

//   template <typename U>
//   friend bool operator==(const Modular<U>& lhs, const Modular<U>& rhs);

//   template <typename U>
//   friend bool operator<(const Modular<U>& lhs, const Modular<U>& rhs);

//   template <typename V, typename U>
//   friend V& operator>>(V& stream, Modular<U>& number);

//  private:
//   Type value;
// };

// template <typename T> bool operator==(const Modular<T>& lhs, const Modular<T>& rhs) { return lhs.value == rhs.value; }
// template <typename T, typename U> bool operator==(const Modular<T>& lhs, U rhs) { return lhs == Modular<T>(rhs); }
// template <typename T, typename U> bool operator==(U lhs, const Modular<T>& rhs) { return Modular<T>(lhs) == rhs; }

// template <typename T> bool operator!=(const Modular<T>& lhs, const Modular<T>& rhs) { return !(lhs == rhs); }
// template <typename T, typename U> bool operator!=(const Modular<T>& lhs, U rhs) { return !(lhs == rhs); }
// template <typename T, typename U> bool operator!=(U lhs, const Modular<T>& rhs) { return !(lhs == rhs); }

// template <typename T> bool operator<(const Modular<T>& lhs, const Modular<T>& rhs) { return lhs.value < rhs.value; }

// template <typename T> Modular<T> operator+(const Modular<T>& lhs, const Modular<T>& rhs) { return Modular<T>(lhs) += rhs; }
// template <typename T, typename U> Modular<T> operator+(const Modular<T>& lhs, U rhs) { return Modular<T>(lhs) += rhs; }
// template <typename T, typename U> Modular<T> operator+(U lhs, const Modular<T>& rhs) { return Modular<T>(lhs) += rhs; }

// template <typename T> Modular<T> operator-(const Modular<T>& lhs, const Modular<T>& rhs) { return Modular<T>(lhs) -= rhs; }
// template <typename T, typename U> Modular<T> operator-(const Modular<T>& lhs, U rhs) { return Modular<T>(lhs) -= rhs; }
// template <typename T, typename U> Modular<T> operator-(U lhs, const Modular<T>& rhs) { return Modular<T>(lhs) -= rhs; }

// template <typename T> Modular<T> operator*(const Modular<T>& lhs, const Modular<T>& rhs) { return Modular<T>(lhs) *= rhs; }
// template <typename T, typename U> Modular<T> operator*(const Modular<T>& lhs, U rhs) { return Modular<T>(lhs) *= rhs; }
// template <typename T, typename U> Modular<T> operator*(U lhs, const Modular<T>& rhs) { return Modular<T>(lhs) *= rhs; }

// template <typename T> Modular<T> operator/(const Modular<T>& lhs, const Modular<T>& rhs) { return Modular<T>(lhs) /= rhs; }
// template <typename T, typename U> Modular<T> operator/(const Modular<T>& lhs, U rhs) { return Modular<T>(lhs) /= rhs; }
// template <typename T, typename U> Modular<T> operator/(U lhs, const Modular<T>& rhs) { return Modular<T>(lhs) /= rhs; }

// template<typename T, typename U>
// Modular<T> power(const Modular<T>& a, const U& b) {
//   assert(b >= 0);
//   Modular<T> x = a, res = 1;
//   U p = b;
//   while (p > 0) {
//     if (p & 1) res *= x;
//     x *= x;
//     p >>= 1;
//   }
//   return res;
// }

// template <typename T>
// bool IsZero(const Modular<T>& number) {
//   return number() == 0;
// }

// template <typename T>
// string to_string(const Modular<T>& number) {
//   return to_string(number());
// }

// // U == std::ostream? but done this way because of fastoutput
// template <typename U, typename T>
// U& operator<<(U& stream, const Modular<T>& number) {
//   return stream << number();
// }

// // U == std::istream? but done this way because of fastinput
// template <typename U, typename T>
// U& operator>>(U& stream, Modular<T>& number) {
//   typename common_type<typename Modular<T>::Type, long long>::type x;
//   stream >> x;
//   number.value = Modular<T>::normalize(x);
//   return stream;
// }

// // using ModType = int;

// // struct VarMod { static ModType value; };
// // ModType VarMod::value;
// // ModType& md = VarMod::value;
// // using Mint = Modular<VarMod>;

// constexpr int md = (int)1e9+7;
// using Mint = Modular<std::integral_constant<decay<decltype(md)>::type, md>>;

// // vector<Mint> fact(1, 1);
// // vector<Mint> inv_fact(1, 1);

// // Mint C(int n, int k) {
// //   if (k < 0 || k > n) {
// //     return 0;
// //   }
// //   while ((int) fact.size() < n + 1) {
// //     fact.push_back(fact.back() * (int) fact.size());
// //     inv_fact.push_back(1 / fact.back());
// //   }
// //   return fact[n] * inv_fact[k] * inv_fact[n - k];
// // }
// int precalculated = 1;
// vector<int> least= {0,1},primes;
// void RunSieve(int n) {
//   n = max(n, 1);
//   least.assign(n + 1, 0);
//  for(int i=4;i<=n;i+=2)least[i]=2;
 
//   for (int i = 3; i * i <= n; i+=2) {
//     if (least[i] == 0) {
//       for (int j = i * i; j <= n; j += i) {
//         if (least[j] == 0) {
//           least[j] = i;
//         }
//       }
//     }
//   }
 
//   for (int i = 2; i <= n; i++) {
//     if (least[i] == 0) {
//       least[i] = i;
//       primes.push_back(i);
//     }
//   }
//   precalculated = n;
// }
// template<typename T> vector<pair<T,int>> Factorize(T n) {
//     if (n <= precalculated) {
//           vector<pair<T, int>> ret;
//           while (n > 1) {
//             if (!ret.empty() && ret.back().first == least[n]) {
//               ret.back().second++;
//             } else {
//               ret.emplace_back(least[n], 1);
//             }
//             n /= least[n];
//           }
//           return ret;
//   }else{ 
//        vector<pair<T,int>> ret;
//           for (int d : {2, 3, 5}) {
//             int cnt=0;
//               while (n % d == 0) {
//                   n /= d;
//                 cnt++;
//               }
//             if(cnt>0)ret.emplace_back(d,cnt); 
//           }
//           static array<int, 8> increments = {4, 2, 4, 2, 4, 6, 2, 6}; int i = 0;
//           for (int64_t d = 7; d * d <= n; d += increments[i++]) {
//             int cnt=0;
//               while (n % d == 0) {
//                   n /= d;cnt++;
//               }
//               if(cnt>0)ret.emplace_back(d,cnt);
//               if(i == 8)i=0;
//           }
//           if (n > 1)
//               ret.emplace_back(n,1);
//           return ret;
//       }
// } 
// template <typename T>int64_t Factorize_to_all(T x) {
// /* num of factors =τ(n)= 1 to k(αi + 1)prduct
// where (αi) is primes power and k is number of primes in prime factorizer;
// sum of factor σ(n) = 1 to k ((pi^(αi+1))-1)/(pi-1) prduct;
// number of co prime for N Euler’s Phi (Totient) function ϕ(N) = N × (1 to k ) (1 − (1/pi); */  
//  if (x == 1)return x;
//   if (x <= precalculated) {
//        int64_t div_sum=1;//for diviors sum
//        int divior_cnt=1;// for count diviors;
//        int coprime_cnt=x;// for coprime cnt;
//        while(x>1){
//            ll factor=least[x], base=least[x];int power=1;
//            while(least[x]==factor){
//                x/=least[x];
//                base*=factor;
//                power++;
//              if(x<=1)break;
//            }
//           coprime_cnt-=(coprime_cnt/factor); 
//           div_sum*=(base-1)/(factor-1);
//           divior_cnt*=power;
//        }
//     return div_sum;
//   }else{
//      assert(false);
//   }
// }
// template <typename T>vector<T> BuildDivisorsFromFactors(const vector<pair<T, int>>& factors) {
//   vector<T> divisors = {1};
//   for (auto& p : factors) {
//     int sz = (int) divisors.size();
//     for (int i = 0; i < sz; i++) {
//       T cur = divisors[i];
//       for (int j = 0; j < p.second; j++) {
//         cur *= p.first;
//         divisors.push_back(cur);
//       }
//     }
//   }
//   sort(divisors.begin(), divisors.end());
//   return divisors;
// }
// vector<ll> Divisors(ll n) {
//  vector<ll> divisors;
//     for (long long i = 1; i * i <= n; ++i) {
//         if (n % i == 0) {
//             divisors.push_back(i);
//             if (n / i != i)
//                 divisors.push_back(n / i);
//         }
//     }
//     return divisors;
// }
// struct Node{
     
// };
// struct segtree{
//   vector<Node>tree;
//    int size;
//    Node NUTURL_ELEMENT={};
//    ll no_oparation ;
//      segtree(int n){
//            assert(n>0);
//             size=1;
//            while(size<n)size<<=1;
//             tree.assign(size<<1,{});
//             size=n-1;
//         }
//         template <class T>     
//         inline void leaf(Node &a,T value){

//         }
//         template<class T>
//         inline void apply(Node &a,T value,int lx=0,int rx=0){
           
//         }
//         inline void push(int x,int lx,int rx){

//         /*
//            if(lx==rx)return;
//            if(tree[x].add!=no_oparation){
//              apply(tree[(x<<1)+1],tree[x].add);
//              apply(tree[(x<<1)+2],tree[x].add);
//              tree[x].add=no_oparation;
//            }

//           */ 
//         }
//        inline Node unite(const Node &a,const  Node &b){
         
          
//         }
//         template<class T>
//        void build (vector<T>&v,int x, int lx,int rx){
//             if(lx==rx && lx<sz(v)){
//            leaf(tree[x],v[lx]);
//             return;
//             }
//         int mid =(lx+rx)>>1;
//             build(v,(x<<1)|1,lx ,mid);
//             build (v,(x<<1)+2,mid+1,rx);

//         tree[x]=unite(tree[(x<<1)|1],tree[(x<<1)+2]);
//         }
//         template<class T>
//         void build(vector<T>&a){
//             build(a,0,0,size);
//         }
//         template<class T>
//         void set(int x , int lx ,int rx, int l,int r, T &value){
//            push(x,lx,rx);
//             if (l > rx || r < lx)
//                 return ; 
//             if (lx >= l && rx <= r){
//                apply(tree[x],value,lx,rx);         
//                 return;
//             }  
//             int mid =(lx+rx)>>1;
//            set((x<<1)|1,lx,mid,l,r,value);
//             set((x<<1)+2,mid+1,rx,l,r,value);
//             tree[x]=unite(tree[(x<<1)|1],tree[(x<<1)+2]);
//         }
//         template<class T>
//         void set(int l ,int r,T &val){
//           assert(0 <= l && l <= r && r <= size);
//             set(0,0,size,l,r,val);
//         }
//         Node get(int x, int lx, int rx, int l, int r){
//             push(x,lx,rx);
//             if (l > rx || r < lx)
//                 return NUTURL_ELEMENT; 
//             if (lx >= l && rx <= r)
//                 return tree[x];

//             int mid =(lx+rx)>>1;
//           Node p1 = get((x<<1)|1,lx,mid,l,r);
//           Node  p2 = get((x<<1)+2,mid+1,rx,l,r);
//             return unite(p1,p2); 
//         }
//        Node get(int l,int r){
//         assert(0 <= l && l <= r && r <= size);
//            return get(0,0,size,l,r);
//         }
        
// };
  // auto good=[&](ll x)->bool {

  //  };
  //  ll l=-1,r=1;
  //   while(!good(r))r<<=2;
    
  //   while(l+1<r){
  //    ll mid=(l+r)>>1;
  //    if(good(mid))r=mid;
  //    else l=mid;
  //   }
  // struct graph{
  //       int n,m,g_time=0;
  //       vector<vector<int>>g;
  //       vector<bool>vis;
  //       vector<int>pa,dis,col,order; 
  //      void resize(){
  //           assert(n>0);
  //           g.assign(n+1,{});
  //           dis.resize(n+1);
  //           pa.resize(n+1,-1);
  //           order.clear();
  //           col.resize(n+1);
  //           vis.resize(n+1);
  //          assert(pa.size()>0 && dis.size()>0 && col.size()>0 && vis.size()>0);
  //      }
  //      void dfs(int s){
  //           order.emplace_back(s);
  //               vis[s]=1;
  //               for(auto i:g[s]){
  //                 if(!vis[i]){
  //                    pa[i]=s;
  //                    dis[i]=dis[s]+1;
  //                    dfs(i);
  //                 }
  //               }
  //      }
  //      void dfs_all(){
  //        for(int i=1;i<=n;i++){ 
  //          if(!vis[i]){
  //            dfs(i);
  //          }
  //        }
  //        assert(order.size()==n);
  //      }
  //      void bfs(int s){
  //             queue<int> q;
  //             q.emplace(s);
  //             vis[s] = true;
  //             pa[s] = -1;
  //              dis[s]=0 ;
  //             while (!q.empty()) {
  //                 int v = q.front();
  //                 q.pop();
  //                 for (int u : g[v]) {
  //                     if (!vis[u]) {
  //                         vis[u] = true;
  //                         q.emplace(u);
  //                         dis[u] = dis[v] + 1;
  //                         pa[u] = v;
  //                     }
  //                 }
  //             }
  //         }
  //   void Test_case(){
    
  //   }
  // };
//  struct edge{
//        int fo,w;
//        edge(int fo,int w):fo(fo),w(w){}
//        bool operator<(edge const &e)const {
//          return e.w<w;
//        }
//  };
//  void  Dijkstra(){
//          int n,m;
//           cin>>n>>m;
//          vector<vector<pair<int,int>>>g(n+1);
//          vector<int>dis(n+1,INT32_MAX),pa(n+1);
         
//           for(int i=0;i<m;i++){ 
//           int u,v,w;cin>>u>>v>>w;
//           g[u].emplace_back(v,w);
//           g[v].emplace_back(u,w);
//          }
//          priority_queue<edge>q;
//          int sur=1;
//          q.emplace(sur,0);
//          dis[sur]=0;
//           while(!q.empty()){
//             auto v=q.top();
//             q.pop();
//             for(auto &[u,w]:g[v.fo]){
//                  if(w+dis[v.fo]<dis[u]){ 
//                    dis[u]=w+dis[v.fo];
//                    q.emplace(u,dis[u]);
//                  }
//             }          
//          }
//  }
// int n,m;
// string girid[1005];
// int level [1005][1005];
// int dx[]={1,-1,0,0};
// int dy []={0,0,1,-1};
// map<pair<int,int> ,pair<int,int>>parant;
// vector<char>Direction;
// bool isvalid(int x, int y){
//      return x>=0 and x<n and y>=0 and y<m and girid[x][y]!='#';
//  }
// void bfs(pair<int ,int> s){
//     queue<pair<int,int>>q;
//      q.push(s);
//      level[s.ff][s.ss]=0;
//      parant[{s.ff,s.ss}]={-1,-1};
//      while(!q.empty()){
//       pii u=q.front();
//       q.pop();
//           for(int i= 0;i<4;i++){
//             int x=u.ff+dx[i];
//             int y=u.ss+dy[i];
//               if(isvalid(x,y) && level[x][y]==-1){
//                   q.push({x,y});
//                   level[x][y]=level[u.ff][u.ss]+1;
//                   parant[{x,y}]={u.ff,u.ss};
//               }
//           }
//      }
//      return ;
// }
// ll Bigmod ( ll a,ll p,ll m ){
//     ll res = 1ll; ll x = a;
//     while (p){
//         if (p&1)res=(res*x) % m;
 
//         x=(x*x)%m;
//       p >>= 1ll;
//     }
//     return res;
// }
// ll gcd(ll a,ll b){
//     if (a == 0ll)    return b;
//     return gcd(b % a, a);
// }
// ll lcm(ll a, ll b){
//     return (a*b) / gcd(a, b);
// } 

#define findc(cn, abc) ((cn).find(abc) != (cn).end())
#define clr(cnt, x) memset((cnt), (x), sizeof(cnt))
template <class T>using V=vector<T>;
#define all(x) (x).begin(), (x).end()
#define lb lower_bound
#define ub upper_bound
#define bg(cn) cn.begin()
template<typename T>
istream &operator>>(istream &in,vector<T>&a){for(T &i:a)in>>i;return in;}
template<typename T>ostream &operator<<(ostream &out, const vector<T> &a) {for (auto i : a) {out << i << " ";}return out;}
template<typename T>void remDup(vector<T> &v){sort(all(v)); v.erase(unique(all(v)),end(v));}
template<typename T,typename U>void safeErase(T &t, const U &u) {auto it = t.find(u);assert(it != end(t)); t.erase(it);}
template<class T>int lwb(const V<T>&a,const T&b){return int(lb(all(a), b)-bg(a));}// (v[r]>=x && v[l]<x) return r;
template<class T>int upb(const V<T>&a,const T&b){return int(ub(all(a),b)-bg(a));}//( v[r]>x && v[l]<=x) return r;
ll SUM_(ll n){return (n*(n+1))/2;}
ll SUM_(ll a,ll b) { return (b- a+1) * (a + b) / 2;}
ll A_SUM(ll n,ll a,ll b){return (n*(a+b))/2;}
ll G_SUM(ll k,ll a,ll b){return ((b*k)-a)/(k-1);}
ll cdiv(ll a, ll b) {return a / b + ((a ^ b) > 0 && a % b);} 
ll fdiv(ll a, ll b) {return a / b - ((a ^ b) < 0 && a % b);} 
const ll inf = (ll)1e18;
const ll _inf=(ll)-1e18;  
const int dx[4]{1, 0, -1, 0},dy[4]{0, 1, 0, -1};  
#define PQ priority_queue
template<class T>using  PQ1 =PQ<T,V<T>,greater<T>>;
template<class T>using  PQ2 =PQ<T,V<pair<T,T>>,greater<pair<T,T>>>; 
#define ff first
#define eb emplace_back
#define ss second 
#define pb push_back
#define to_str(x) (to_string(x)) 
#define to_ll(x) (stoi(x))
#define to_char(x) ((x)+'0')
#define to_int(x) ((x)-'0')
constexpr int Pct(int x){return __builtin_popcount(x);}  
constexpr int Msb(int x){return x==0?0:31 -__builtin_clz(x);}  
constexpr int P2(int x) { return 1 << x; }
constexpr int Msk2(int x) { return P2(x) - 1;}
constexpr int Lsb(int x) { return  x&~(x-1);}
typedef vector<pair<ll,ll>> vll;
typedef vector<pair<int,int>>vii;
typedef vector<ll>vl;  
typedef pair<int,int> pii;
typedef pair<ll, ll> pll;
#define sz(n) (int)(n).size() 
#define fora(cn) for(auto &x : (cn))
#define PI 2*acos(0) 
#define el '\n'
typedef unsigned long long ull;
typedef vector<int> vi;
#define se(X) setprecision(X)
#define rall(x) (x).rbegin(), (x).rend()
#define rep(n) for(int (i)=(0);(i)<(n);(i)++)
#define forr(z,x,n) for(int (z)=(x);(z)<(n);(z)++)
#define rfor(z,x,n) for(int (z)=(x);(z)>=(n);(z)--)
const int MX = (int)2e5 + 5;
int main(){ 
  ios::sync_with_stdio(false);cin.tie(0);
     
}
