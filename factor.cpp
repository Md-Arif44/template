
int precalculated = 1;
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
  precalculated = n;
}
template<typename T> vector<pair<T,int>> Factorize(T n) {
    if (n <= precalculated) {
          vector<pair<T, int>> ret;
          while (n > 1) {
            if (!ret.empty() && ret.back().first == least[n]) {
              ret.back().second++;
            } else {
              ret.emplace_back(least[n], 1);
            }
            n /= least[n];
          }
          return ret;
  }else{ 
       vector<pair<T,int>> ret;
          for (int d : {2, 3, 5}) {
            int cnt=0;
              while (n % d == 0) {
                  n /= d; cnt++;
              }
            if(cnt>0)ret.emplace_back(d,cnt); 
          }
          static array<int, 8> increments = {4, 2, 4, 2, 4, 6, 2, 6}; int i = 0;
          for (int64_t d = 7; d * d <= n; d += increments[i++]) {
            int cnt=0;
              while (n % d == 0) {
                  n /= d;cnt++;
              }
              if(cnt>0)ret.emplace_back(d,cnt);
              if(i == 8)i=0;
          }
          if (n > 1)
              ret.emplace_back(n,1);
          return ret;
      }
} 
template <typename T>int64_t Factorize_to_all(T x) {
/* num of factors =τ(n)= 1 to k(αi + 1)prduct
where (αi) is primes power and k is number of primes in prime factorizer;
sum of factor σ(n) = 1 to k ((pi^(αi+1))-1)/(pi-1) prduct;
number of co prime for N Euler’s Phi (Totient) function ϕ(N) = N × (1 to k ) (1 − (1/pi); */  
 if (x == 1)return x;
  if (x <= precalculated) {
       int64_t div_sum=1;//for diviors sum
       int divior_cnt=1;// for count diviors;
       int coprime_cnt=x;// for coprime cnt;
       while(x>1){
           ll factor=least[x], base=least[x];int power=1;
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
template <typename T>vector<T> BuildDivisorsFromFactors(const vector<pair<T, int>>& factors) {
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
