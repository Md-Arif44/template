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
  T sum(int l,int r ) { // [l,r) 
    return sum(r) - sum(l) ;}
};
