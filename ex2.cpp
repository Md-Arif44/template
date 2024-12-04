ll Bigmod ( ll a,ll p,ll m ){
    ll res = 1ll; ll x = a;
    while (p){
        if (p&1)res=(res*x) % m;
 
        x=(x*x)%m;
      p >>= 1ll;
    }
    return res;
}
// div
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
ll gcd(ll a,ll b){
    if (a == 0ll)    return b;
    return gcd(b % a, a);
}
ll lcm(ll a, ll b){
    return (a*b) / gcd(a, b);
} 
//bs
  auto good=[](ll x)->bool {

   };
   ll l=-1,r=1;
    while(!good(r))r<<=2;
    
    while(l+1<r){
     ll mid=(l+r)>>1;
     if(good(mid))r=mid;
     else l=mid;
    }
//Queue  
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
