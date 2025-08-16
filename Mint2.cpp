template< const int  m >
 struct Mint{
  int x;
  Mint() : x(0){}
  Mint(long long x_):x(x_ % m){if (x < 0) x += m;}
  template <typename U>explicit operator U() const { return static_cast<U>(x); }
  Mint &operator+=(Mint b){if ((x += b.x) >= m) x -= m; return *this;}
  Mint &operator-=(Mint b){if ((x -= b.x) < 0) x += m; return *this;}
  Mint &operator*=(Mint b){x= (long long)(x) * b.x % m; return *this;}
  Mint power(long long e) const {
      Mint r = 1,b =*this;
     for (;e;e>>=1,b*=b) if(e&1)r*=b;
      return r;
  }
  Mint inv(){return power(m - 2);}
  Mint &operator/=(Mint b){return *this *= b.power(m - 2);}
  friend Mint operator+(Mint a, Mint b){return a += b;}
  friend Mint operator-(Mint a, Mint b){return a -= b;}
  friend Mint operator/(Mint a, Mint b){return a /= b;}
  friend Mint operator*(Mint a, Mint b){return a *= b;}
  friend bool operator==(Mint a, Mint b){return a.x == b.x;}
  friend bool operator!=(Mint a, Mint b){return a.x != b.x;}
  friend std::ostream &operator<<(std::ostream &os, Mint const &a) { return os << a.x; }
  friend std::istream& operator>>(std::istream &is, Mint& a) {is >> a.x;a.x%= m;if(a.x< 0)a.x+=m;return is;}
 };
 // using mint = Mint< 998244353 > ;
 using mint = Mint< int(1e9)+7 > ;
 
 vector<mint> fact(1, 1);
 vector<mint> inv_fact(1, 1);
 
 mint C(int n, int k) {
  if (k < 0 || k > n) return 0;
  
  while ((int) fact.size() < n + 1) {
    fact.push_back(fact.back() * (int) fact.size());
    inv_fact.push_back( mint(1) / fact.back());
  }
  return fact[n] * inv_fact[k] * inv_fact[n - k];

}
