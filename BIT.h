template <class T=ll >
struct FT {
  vector<T> x;
  FT(int n) : x(n+1) { } 
  void add(int k, T a) { // x[k] += a
    for (++k; k < x.size(); k += k&-k) x[k] += a;
  }
  T sum(int k,T s=0) { // [0,k)
    for (; k > 0; k &= k-1) s += x[k];
    return s;
  }
  T rsum(int l,int r ) { // [l,r) 
    return sum(r) - sum(l) ;}
};


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
