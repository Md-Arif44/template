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
// subarray
  vector<int> lop(n),rop(n) ,stk;  
    rep(i,n){
        while(sz(stk) && a[stk.back()]<a[i] ){
             stk.pop_back() ;
        }
        lop[i]= (stk.empty()? -1: stk.back()) ;
       stk.pb(i) ;
    } 
    stk.clear();
    repr(i,n-1,0){
        while(sz(stk) &&  a[stk.back()] <=a[i] ){
             stk.pop_back() ;
        }
      rop[i]= (stk.empty()?n:stk.back() );
      stk.pb(i);
    }

//bs minium  true r
  auto good=[](ll x)->bool {

   };
   ll l=-1,r=1;
    while(!good(r))r<<=2;
    
    while(r >l+1){
     ll mid=(l+r)>>1;
     if(good(mid))r=mid;
     else l=mid;
    }


// bs maximum true l
  auto good=[](ll x)->bool {

   };
   ll l=0,r=1;
    while(good(r))r<<=2;
    
    while(r >l+1){
     ll mid=(l+r)>>1;
     if(good(mid))l=mid;
     else r=mid;
    }


    
vvl PrefixSums2D(vvl &a) {
  // ...1 base
    ll n = a.size(),m = a[0].size();
    vvl prefixSum(n + 1, vl(m + 1, 0));
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
