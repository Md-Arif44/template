// bs
  auto good=[](ll x)->bool {

   };
   ll l=-1,r=1;
    while(!good(r))r<<=2;
    
    while(l+1<r){
     ll mid=(l+r)>>1;
     if(good(mid))r=mid;
     else l=mid;
    }
//presum
template<class T>
V<T> pre_sum(vector<T>&a){
  int X=sz(a);V<T>pre(X);
    for(int i=0;i<X;i++){
      ((i==0)?pre[i]=a[i]:pre[i]=pre[i-1]+a[i]);    
    }
  return pre;
}
template<class T>
V<T> suf_sum(vector<T>&a){
  int X=sz(a);V<T> suf(X);
        for(int i=X-1;i>=0;i--){
          ((i==X-1)?suf[i]=a[i]:suf[i]=suf[i+1]+a[i]);
        }
  return suf;
}
template<class T>
T get(int l,int r,vector<T>&pre){
  assert(0<=l && l<=r &&r<=pre.size()-1);
  if(l==0)return pre[r];
  else return pre[r]-pre[l-1];
}
//queue
struct Two_stack{
        vector<ll>s,op1,op2;
       void push(ll x){
         s.push_back(x);
       }           
       ll remove(){
         ll res=s.back();
          s.pop_back();
         return res;
       }
       pair<ll,ll> get(){
        return {};
       }
       int size(){
        return s.size();
       }
};
Two_stack pre,suf;
void add(ll x){
   suf.push(x);
}
void remove(){
    if(!pre.size()){
         while(suf.size()){
             ll x=suf.remove();
             pre.push(x);
         }
    }
    pre.remove();
}
bool good(){
  
}
// diviors
vector<ll> Divisors(ll n) {
 vector<ll> divisors;
    for (long long i = 1; i * i <= n; ++i) {
        if (n % i == 0) {
            divisors.push_back(i);
            if (n / i != i)
                divisors.push_back(n / i);
        }
    }
    return divisors;
}
// Mint 
template <typename T>
T inverse(T a, T m) {
  T u = 0, v = 1;
  while (a != 0) {
    T t = m / a;
    m -= t * a; swap(a, m);
    u -= t * v; swap(u, v);
  }
  assert(m == 1);
  return u;
}

template <typename T>
class Modular {
 public:
  using Type = typename decay<decltype(T::value)>::type;

  constexpr Modular() : value() {}
  template <typename U>
  Modular(const U& x) {
    value = normalize(x);
  }

  template <typename U>
  static Type normalize(const U& x) {
    Type v;
    if (-mod() <= x && x < mod()) v = static_cast<Type>(x);
    else v = static_cast<Type>(x % mod());
    if (v < 0) v += mod();
    return v;
  }

  const Type& operator()() const { return value; }
  template <typename U>
  explicit operator U() const { return static_cast<U>(value); }
  constexpr static Type mod() { return T::value; }

  Modular& operator+=(const Modular& other) { if ((value += other.value) >= mod()) value -= mod(); return *this; }
  Modular& operator-=(const Modular& other) { if ((value -= other.value) < 0) value += mod(); return *this; }
  template <typename U> Modular& operator+=(const U& other) { return *this += Modular(other); }
  template <typename U> Modular& operator-=(const U& other) { return *this -= Modular(other); }
  Modular& operator++() { return *this += 1; }
  Modular& operator--() { return *this -= 1; }
  Modular operator++(int) { Modular result(*this); *this += 1; return result; }
  Modular operator--(int) { Modular result(*this); *this -= 1; return result; }
  Modular operator-() const { return Modular(-value); }

  template <typename U = T>
  typename enable_if<is_same<typename Modular<U>::Type, int>::value, Modular>::type& operator*=(const Modular& rhs) {
    value = normalize(static_cast<int64_t>(value) * static_cast<int64_t>(rhs.value));
    return *this;
  }
  template <typename U = T>
  typename enable_if<is_same<typename Modular<U>::Type, long long>::value, Modular>::type& operator*=(const Modular& rhs) {
    long long q = static_cast<long long>(static_cast<long double>(value) * rhs.value / mod());
    value = normalize(value * rhs.value - q * mod());
    return *this;
  }
  template <typename U = T>
  typename enable_if<!is_integral<typename Modular<U>::Type>::value, Modular>::type& operator*=(const Modular& rhs) {
    value = normalize(value * rhs.value);
    return *this;
  }

  Modular& operator/=(const Modular& other) { return *this *= Modular(inverse(other.value, mod())); }

  friend const Type& abs(const Modular& x) { return x.value; }

  template <typename U>
  friend bool operator==(const Modular<U>& lhs, const Modular<U>& rhs);

  template <typename U>
  friend bool operator<(const Modular<U>& lhs, const Modular<U>& rhs);

  template <typename V, typename U>
  friend V& operator>>(V& stream, Modular<U>& number);

 private:
  Type value;
};

template <typename T> bool operator==(const Modular<T>& lhs, const Modular<T>& rhs) { return lhs.value == rhs.value; }
template <typename T, typename U> bool operator==(const Modular<T>& lhs, U rhs) { return lhs == Modular<T>(rhs); }
template <typename T, typename U> bool operator==(U lhs, const Modular<T>& rhs) { return Modular<T>(lhs) == rhs; }

template <typename T> bool operator!=(const Modular<T>& lhs, const Modular<T>& rhs) { return !(lhs == rhs); }
template <typename T, typename U> bool operator!=(const Modular<T>& lhs, U rhs) { return !(lhs == rhs); }
template <typename T, typename U> bool operator!=(U lhs, const Modular<T>& rhs) { return !(lhs == rhs); }

template <typename T> bool operator<(const Modular<T>& lhs, const Modular<T>& rhs) { return lhs.value < rhs.value; }

template <typename T> Modular<T> operator+(const Modular<T>& lhs, const Modular<T>& rhs) { return Modular<T>(lhs) += rhs; }
template <typename T, typename U> Modular<T> operator+(const Modular<T>& lhs, U rhs) { return Modular<T>(lhs) += rhs; }
template <typename T, typename U> Modular<T> operator+(U lhs, const Modular<T>& rhs) { return Modular<T>(lhs) += rhs; }

template <typename T> Modular<T> operator-(const Modular<T>& lhs, const Modular<T>& rhs) { return Modular<T>(lhs) -= rhs; }
template <typename T, typename U> Modular<T> operator-(const Modular<T>& lhs, U rhs) { return Modular<T>(lhs) -= rhs; }
template <typename T, typename U> Modular<T> operator-(U lhs, const Modular<T>& rhs) { return Modular<T>(lhs) -= rhs; }

template <typename T> Modular<T> operator*(const Modular<T>& lhs, const Modular<T>& rhs) { return Modular<T>(lhs) *= rhs; }
template <typename T, typename U> Modular<T> operator*(const Modular<T>& lhs, U rhs) { return Modular<T>(lhs) *= rhs; }
template <typename T, typename U> Modular<T> operator*(U lhs, const Modular<T>& rhs) { return Modular<T>(lhs) *= rhs; }

template <typename T> Modular<T> operator/(const Modular<T>& lhs, const Modular<T>& rhs) { return Modular<T>(lhs) /= rhs; }
template <typename T, typename U> Modular<T> operator/(const Modular<T>& lhs, U rhs) { return Modular<T>(lhs) /= rhs; }
template <typename T, typename U> Modular<T> operator/(U lhs, const Modular<T>& rhs) { return Modular<T>(lhs) /= rhs; }

template<typename T, typename U>
Modular<T> power(const Modular<T>& a, const U& b) {
  assert(b >= 0);
  Modular<T> x = a, res = 1;
  U p = b;
  while (p > 0) {
    if (p & 1) res *= x;
    x *= x;
    p >>= 1;
  }
  return res;
}

template <typename T>
bool IsZero(const Modular<T>& number) {
  return number() == 0;
}

template <typename T>
string to_string(const Modular<T>& number) {
  return to_string(number());
}

// U == std::ostream? but done this way because of fastoutput
template <typename U, typename T>
U& operator<<(U& stream, const Modular<T>& number) {
  return stream << number();
}

// U == std::istream? but done this way because of fastinput
template <typename U, typename T>
U& operator>>(U& stream, Modular<T>& number) {
  typename common_type<typename Modular<T>::Type, long long>::type x;
  stream >> x;
  number.value = Modular<T>::normalize(x);
  return stream;
}

// using ModType = int;

// struct VarMod { static ModType value; };
// ModType VarMod::value;
// ModType& md = VarMod::value;
// using Mint = Modular<VarMod>;

constexpr int md = (int)1e9+7;
using Mint = Modular<std::integral_constant<decay<decltype(md)>::type, md>>;

// vector<Mint> fact(1, 1);
// vector<Mint> inv_fact(1, 1);

// Mint C(int n, int k) {
//   if (k < 0 || k > n) {
//     return 0;
//   }
//   while ((int) fact.size() < n + 1) {
//     fact.push_back(fact.back() * (int) fact.size());
//     inv_fact.push_back(1 / fact.back());
//   }
//   return fact[n] * inv_fact[k] * inv_fact[n - k];
// }


//factorize
namespace factorizer {

template <typename T>
struct FactorizerVarMod { static T value; };
template <typename T>
T FactorizerVarMod<T>::value;

template <typename T>
bool IsPrime(T n, const vector<T>& bases) {
  if (n < 2) {
    return false;
  }
  vector<T> small_primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29};
  for (const T& x : small_primes) {
    if (n % x == 0) {
      return n == x;
    }
  }
  if (n < 31 * 31) {
    return true;
  }
  int s = 0;
  T d = n - 1;
  while ((d & 1) == 0) {
    d >>= 1;
    s++;
  }
  FactorizerVarMod<T>::value = n;
  for (const T& a : bases) {
    if (a % n == 0) {
      continue;
    }
    Modular<FactorizerVarMod<T>> cur = a;
    cur = power(cur, d);
    if (cur == 1) {
      continue;
    }
    bool witness = true;
    for (int r = 0; r < s; r++) {
      if (cur == n - 1) {
        witness = false;
        break;
      }
      cur *= cur;
    }
    if (witness) {
      return false;
    }
  }
  return true;
}

bool IsPrime(int64_t n) {
  return IsPrime(n, {2, 325, 9375, 28178, 450775, 9780504, 1795265022});
}

bool IsPrime(int32_t n) {
  return IsPrime(n, {2, 7, 61});
}

// but if you really need uint64_t version...
/*
bool IsPrime(uint64_t n) {
  if (n < 2) {
    return false;
  }
  vector<uint32_t> small_primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29};
  for (uint32_t x : small_primes) {
    if (n == x) {
      return true;
    }
    if (n % x == 0) {
      return false;
    }
  }
  if (n < 31 * 31) {
    return true;
  }
  uint32_t s = __builtin_ctzll(n - 1);
  uint64_t d = (n - 1) >> s;
  function<bool(uint64_t)> witness = [&n, &s, &d](uint64_t a) {
    uint64_t cur = 1, p = d;
    while (p > 0) {
      if (p & 1) {
        cur = (__uint128_t) cur * a % n;
      }
      a = (__uint128_t) a * a % n;
      p >>= 1;
    }
    if (cur == 1) {
      return false;
    }
    for (uint32_t r = 0; r < s; r++) {
      if (cur == n - 1) {
        return false;
      }
      cur = (__uint128_t) cur * cur % n;
    }
    return true;
  };
  vector<uint64_t> bases_64bit = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
  for (uint64_t a : bases_64bit) {
    if (a % n == 0) {
      return true;
    }
    if (witness(a)) {
      return false;
    }
  }
  return true;
}
*/

vector<int> least = {0, 1};
vector<int> primes;
int precalculated = 1;

void RunLinearSieve(int n) {
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
  precalculated = n;
}

void RunSlowSieve(int n) {
  
  n = max(n, 1);
  least.assign(n + 1, 0);
  
 for(int i=4;i<=n;i+=2)least[i]=2;
 
  for (int i = 3; i * i <= n; i+=2) {
    if (least[i] == 0) {
      for (int j = i * i; j <= n; j += i) {
        if (least[j] == 0) {
          least[j] = i;
        }
      }
    }
  }
 
  for (int i = 2; i <= n; i++) {
    if (least[i] == 0) {
      least[i] = i;
      primes.push_back(i);
    }
  }
  precalculated = n;
}

void RunSieve(int n) {
  RunLinearSieve(n);
}

template <typename T>
vector<pair<T, int>> MergeFactors(const vector<pair<T, int>>& a, const vector<pair<T, int>>& b) {
  vector<pair<T, int>> c;
  int i = 0;
  int j = 0;
  while (i < (int) a.size() || j < (int) b.size()) {
    if (i < (int) a.size() && j < (int) b.size() && a[i].first == b[j].first) {
      c.emplace_back(a[i].first, a[i].second + b[j].second);
      ++i;
      ++j;
      continue;
    }
    if (j == (int) b.size() || (i < (int) a.size() && a[i].first < b[j].first)) {
      c.push_back(a[i++]);
    } else {
      c.push_back(b[j++]);
    }
  }
  return c;
}

template <typename T>
vector<pair<T, int>> RhoC(const T& n, const T& c) {
  if (n <= 1) {
    return {};
  }
  if ((n & 1) == 0) {
    return MergeFactors({{2, 1}}, RhoC(n / 2, c));
  }
  if (IsPrime(n)) {
    return {{n, 1}};
  }
  FactorizerVarMod<T>::value = n;
  Modular<FactorizerVarMod<T>> x = 2;
  Modular<FactorizerVarMod<T>> saved = 2;
  T power = 1;
  T lam = 1;
  while (true) {
    x = x * x + c;
    T g = __gcd((x - saved)(), n);
    if (g != 1) {
      return MergeFactors(RhoC(g, c + 1), RhoC(n / g, c + 1));
    }
    if (power == lam) {
      saved = x;
      power <<= 1;
      lam = 0;
    }
    lam++;
  }
  return {};
}

template <typename T>
vector<pair<T, int>> Rho(const T& n) {
  return RhoC(n, static_cast<T>(1));
}

template <typename T>
vector<pair<T, int>> Factorize(T x) {
  if (x <= 1) {
    return {};
  }
  if (x <= precalculated) {
    vector<pair<T, int>> ret;
    while (x > 1) {
      if (!ret.empty() && ret.back().first == least[x]) {
        ret.back().second++;
      } else {
        ret.emplace_back(least[x], 1);
      }
      x /= least[x];
    }
    return ret;
  }
  if (x <= static_cast<int64_t>(precalculated) * precalculated) {
    vector<pair<T, int>> ret;
    if (!IsPrime(x)) {
      for (T i : primes) {
        T t = x / i;
        if (i > t) {
          break;
        }
        if (x == t * i) {
          int cnt = 0;
          while (x % i == 0) {
            x /= i;
            cnt++;
          }
          ret.emplace_back(i, cnt);
          if (IsPrime(x)) {
            break;
          }
        }
      }
    }
    if (x > 1) {
      ret.emplace_back(x, 1);
    }
    return ret;
  }
  return Rho(x);
}
template <typename T>
ll Factorize_to_all(T x) {
/* num of factors =τ(n)= 1 to k(αi + 1)prduct
where (αi) is primes power and k is number of primes in prime factorizer;
sum of factor σ(n) = 1 to k ((pi^(αi+1))-1)/(pi-1) prduct;
number of co prime for N;
Euler’s Phi (Totient) function ϕ(N) = N × (1 to k ) (1 − (1/pi);
*/  
 if (x == 1) {
    return x;
  }
  if (x <= precalculated) {
       ll div_sum=1;//for diviors sum
       int divior_cnt=1;// for count diviors;
       int coprime_cnt=x;// for coprime cnt;
       while(x>1){
           ll factor=least[x];
           ll base=least[x];int power=1;
           while(least[x]==factor){
               x/=least[x];
               base*=factor;
               power++;
             if(x<=1)break;
           }
          coprime_cnt-=(coprime_cnt/factor); 
          div_sum*=(base-1)/(factor-1);
          divior_cnt*=power;
       }
    return div_sum;
  }else{
     assert(false);
  }
}
template <typename T>
vector<T> BuildDivisorsFromFactors(const vector<pair<T, int>>& factors) {
  vector<T> divisors = {1};
  for (auto& p : factors) {
    int sz = (int) divisors.size();
    for (int i = 0; i < sz; i++) {
      T cur = divisors[i];
      for (int j = 0; j < p.second; j++) {
        cur *= p.first;
        divisors.push_back(cur);
      }
    }
  }
  sort(divisors.begin(), divisors.end());
  return divisors;
}

}  
using namespace factorizer; 

// lde

template<typename T>
T extgcd(T a, T b, T &x, T &y)
{
    if (a == 0) return x = 0, y = 1, b;
    T p = b / a;
    T g = extgcd(b - p * a, a, y, x);
    x -= p * y;
    return g;
}
/// Return true if there exist such (x, y) satisfy ax + by = c
/// Find (&g) = gcd(a, b)
/// Find (&x, &y) satisfy ax + by = c
template<typename T>
bool find_any_solution(T a, T b, T c, T &x, T &y, T &g)
{
    if (a == 0 && b == 0) /// 0x + 0y = c
    {
        if (c != 0) return false;
        x = y = g = 0;
        return true;
    }

    if (a == 0) /// 0x + by = c
    {
        if (c % b != 0) return false;
        x = 0, y = c / b, g = abs(b);
        return true;
    }

    if (b == 0) /// ax + 0y = c
    {
        if (c % a != 0) return false;
        x = c / a, y = 0, g = abs(a);
        return true;
    }

    /// ax + by = c
    g = extgcd(abs(a), abs(b), x, y);
    if (c % g != 0) return false;

    x *= (a < 0 ? -1 : +1) * c / g;
    y *= (b < 0 ? -1 : +1) * c / g;
    return true;
}
/// Find the next/prev (cnt)-th solution of ax + by = c
template<typename T>
void shift_solution(T & x, T & y, T a, T b, T cnt)
{
    x += cnt * b;
    y -= cnt * a;
}
template<typename T = long long>
T find_all_solutions(T a, T b, T c, T min_x, T max_x, T min_y, T max_y) {
    if (min_x > max_x) return 0; /// Invalid range
    if (min_y > max_y) return 0; /// Invalid range

    if (a == 0 && b == 0) /// 0x + 0y = c
    {
        if (c != 0) return 0; /// No solution
        return 1LL * (max_x - min_x + 1) * (max_y - min_y + 1); /// Ways to select (x) and (y) in range
    }

    if (a == 0) /// 0x + by = c <=> y = c / b
    {
        if (c % b != 0) return 0; /// No solution
        if (1LL * min_y * b > c) return 0; /// Out of range: min > y
        if (1LL * max_y * b < c) return 0; /// Out of range: max < y
        return max_x - min_x + 1; /// Ways to select (x) in range    
    }

    if (b == 0) /// ax + 0y = c <=> x = c / a
    {
        if (c % a != 0) return 0; /// No solution
        if (1LL * min_x * a > c) return 0; /// Out of range: min > x
        if (1LL * max_x * a < c) return 0; /// Out of range: max < x
        return max_y - min_y + 1; /// Ways to select (y) in range    
    }

    T x, y, g;
    if (!find_any_solution(a, b, c, x, y, g)) return 0;
    a /= g;     
    b /= g;

    T sign_a = a > 0 ? +1 : -1;
    T sign_b = b > 0 ? +1 : -1;

    shift_solution(x, y, a, b, (min_x - x) / b);
    if (x < min_x) shift_solution(x, y, a, b, sign_b);
    if (x > max_x) return 0;
    T lx1 = x;

    shift_solution(x, y, a, b, (max_x - x) / b);
    if (x > max_x) shift_solution(x, y, a, b, -sign_b);
    T rx1 = x;

    shift_solution(x, y, a, b, -(min_y - y) / a);
    if (y < min_y) shift_solution(x, y, a, b, -sign_a);
    if (y > max_y) return 0;
    T lx2 = x;

    shift_solution(x, y, a, b, -(max_y - y) / a);
    if (y > max_y) shift_solution(x, y, a, b, sign_a);
    T rx2 = x;

    if (lx2 > rx2)swap(lx2, rx2);
    T lx = max(lx1, lx2);
    T rx = min(rx1, rx2);

    if (lx > rx) return 0;
    return (rx - lx) / abs(b) + 1;
}
Find the solution with minimum value x+y 
 x^=x+k (b/g);
  y^=y-k (a/g);
  if a<b select smalest possible k
  if b>a	select largest possible k
