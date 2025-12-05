
vi suffix_array(int n, const string&s, int char_bound= 260) {
  vi a(n);
  if (n == 0)return a;
  if (char_bound != -1 ) {
    vi aux(char_bound, 0); int sum = 0;
    rep(i,0,n)aux[s[i]]++;
    rep(i,0,char_bound) {
      int add = aux[i];
      aux[i] = sum; sum += add;
    }
    rep( i,0,n) a[aux[s[i]]++] = i;
  } else {
    iota(all(a), 0);
    sort(all(a), [&s](int i, int j) { return s[i] < s[j]; });
  }
  vi sbs(n),  ptg(n),neg(n), gp(n);
  gp[a[0]] = 0;
  rep (i,1,n) {
    gp[a[i]] = gp[a[i - 1]] + (!(s[a[i]] == s[a[i - 1]]));
  }
  int cnt = gp[a[n - 1]] + 1 , step = 1;
  while (cnt < n) {
    int at = 0;
    rep (i,n - step,n )sbs[at++] = i;
    
    rep(i,0,n) {
      if (a[i] - step >= 0) {
        sbs[at++] = a[i] - step;
      }
    }
    for (int i = n - 1; i >= 0; i--)ptg[gp[a[i]]] = i;
    rep(i ,0,n) {
      int x = sbs[i];
      a[ptg[gp[x]]++] = x;
    }
    neg [a[0]] = 0;
    rep(i,1,n) {
      if (gp[a[i]] != gp[a[i - 1]]) {
        neg [a[i]] = neg[a[i - 1]] + 1;
      } else {
        int pre = (a[i - 1] + step >= n ? -1 : gp[a[i - 1] + step]);
        int cur = (a[i] + step >= n ? -1 : gp[a[i] + step]);
        neg[a[i]] = neg[a[i - 1]] + (pre != cur);
      }
    }
    swap(gp, neg);
    cnt = gp[a[n - 1]] + 1;step <<= 1;
  } return a;
}

vi build_lcp(int n, const string &s, const vi &sa) {
  assert((int) sa.size() == n);
  vi pos(n), lcp(max(n - 1, 0)) ;
  rep (i,0,n)pos[sa[i]] = i; int k=0 ;
  rep(i,0,n) {
    k = max(k - 1, 0);
    if (pos[i] == n - 1)k = 0;
    else {
      int j = sa[pos[i] + 1];
      while (i + k < n && j + k < n && s[i + k] == s[j + k]) {
        k++;
      }
      lcp[pos[i]] = k;
    }
  } return lcp;
}
 
