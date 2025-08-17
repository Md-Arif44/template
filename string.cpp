 

struct  Hash{
      int n;
      const ll p1=int(1e9+7),p2=int(1e9+9),b1=int(1e7+7) ,b2=int(1e7+11);
      vector<ll> pw1,pw2,H1,H2 ;
      template<typename T>
      Hash(const T &s){
          n=int(s.size());
          H1.resize(n);H2.resize(n);
          H1[0]=s[0];H2[0]=s[0];
          pw1.resize(n,1);pw2.resize(n,1) ;
          for(int i=1;i<n;i++){
              H1[i]= (H1[i-1]*b1 +s[i]) %p1  ;
              H2[i]= (H2[i-1]*b2 +s[i]) %p2  ;
              pw1[i]= pw1[i-1]*b1% p1 ;     
              pw2[i]= pw2[i-1]*b2% p2 ;      
          }
      }
     inline pair<ll,ll> operator()(const int l,const int r)const {
         assert(l>=0 && l<=r && r<=n-1);  
          if(l==0)return {H1[r],H2[r]};
          else { 
            ll val1= H1[r]-H1[l-1]*pw1[r-l+1] % p1 ; if(val1<0)val1+=p1 ;
            ll val2= H2[r]-H2[l-1]*pw2[r-l+1] % p2 ; if(val2<0)val2+=p2 ;
          return {val1,val2};
         }
      }
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
