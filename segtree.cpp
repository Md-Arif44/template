template <class T, T (*op)(T, T), T (*e)()> 
struct segtree {
  int n,size;
  vector<T> d;
  segtree(const vector<T>& v){
      n=int(v.size()),size=1; 
      while(size<n)size<<=1;
      
        d=vector<T>(2 * size, e());
        for(int i=0;i<n;i++)d[size + i] = v[i];
        for(int i=size-1;i>=1;i--) update(i);
    }
    void update(int k) { 
      d[k] = op(d[2 * k], d[2 * k + 1]); 
    }
    void set(int p, T x) {
        assert(0 <= p && p < n);
        p += size;  d[p] = x;
        while(p > 1)update(p >>=1);
    }
    T prod(int l, int r) const { // [l,r)
      assert(0 <= l && l <= r && r <= n);
        T sml = e(), smr = e();
        for(l += size ,r += size ; l<r ;  l>>=1,r>>=1 ){
            if (l & 1) sml = op(sml, d[l++]);
            if (r & 1) smr = op(d[--r], smr); 
        }
      return op(sml, smr);
    }
};

ll op(ll x ,ll y){
   return x+y;
}
ll e(){
   return 0ll ;
}




template <class S,
          S (*op)(S, S),
          S (*e)(),
          class F,
          S (*Map)(F, S),
          F (*Com)(F, F),
          F (*id)()>
struct lazy_segtree {

    lazy_segtree(int n) : lazy_segtree(std::vector<S>(n, e())) {}
    lazy_segtree(const std::vector<S>& v) : _n(int(v.size())) {  
        size = 1 ;while(size<_n)size*=2 ;
        log = 31 - __builtin_clz(size); 
        d = std::vector<S>(2 * size, e());
        lz = std::vector<F>(size, id());
        for (int i = 0; i < _n; i++) d[size + i] = v[i];
        for (int i = size - 1; i >= 1; i--) update(i);
        
    }
 
    void set(int p, S x) {
        assert(0 <= p && p < _n);
        p += size;
        for (int i = log; i >= 1; i--) push(p >> i);
        d[p] = x;
        for (int i = 1; i <= log; i++) update(p >> i);
    }
 
    S get(int p) {
        assert(0 <= p && p < _n);
        p += size;
        for (int i = log; i >= 1; i--) push(p >> i);
        return d[p];
    }
 
    S prod(int l, int r) { // [l,r)
        assert(0 <= l && l <= r && r <= _n);
        if (l == r) return e();
 
        l += size;
        r += size;
 
        for (int i = log; i >= 1; i--) {
            if (((l >> i) << i) != l) push(l >> i);
            if (((r >> i) << i) != r) push((r - 1) >> i);
        }
 
        S sml = e(), smr = e();
        while (l < r) {
            if (l & 1) sml = op(sml, d[l++]);
            if (r & 1) smr = op(d[--r], smr);
            l >>= 1;
            r >>= 1;
        }
 
        return op(sml, smr);
    }
 
    S all_prod() { return d[1]; }
 
    void apply(int p, F f) {
        assert(0 <= p && p < _n);
        p += size;
        for (int i = log; i >= 1; i--) push(p >> i);
        d[p] = Map(f, d[p]);
        for (int i = 1; i <= log; i++) update(p >> i);
    }
    void apply(int l, int r, F f) { // [l,r)
        assert(0 <= l && l <= r && r <= _n);
        if (l == r) return;
 
        l += size;
        r += size;
 
        for (int i = log; i >= 1; i--) {
            if (((l >> i) << i) != l) push(l >> i);
            if (((r >> i) << i) != r) push((r - 1) >> i);
        }
 
        {
            int l2 = l, r2 = r;
            while (l < r) {
                if (l & 1) all_apply(l++, f);
                if (r & 1) all_apply(--r, f);
                l >>= 1;
                r >>= 1;
            }
            l = l2;
            r = r2;
        }
 
        for (int i = 1; i <= log; i++) {
            if (((l >> i) << i) != l) update(l >> i);
            if (((r >> i) << i) != r) update((r - 1) >> i);
        }
    }
 
    template <bool (*g)(S)> int max_right(int l) {
        return max_right(l, [](S x) { return g(x); });
    }
    template <class G> int max_right(int l, G g) {
        assert(0 <= l && l <= _n);
        assert(g(e()));
        if (l == _n) return _n;
        l += size;
        for (int i = log; i >= 1; i--) push(l >> i);
        S sm = e();
        do {
            while (l % 2 == 0) l >>= 1;
            if (!g(op(sm, d[l]))) {
                while (l < size) {
                    push(l);
                    l = (2 * l);
                    if (g(op(sm, d[l]))) {
                        sm = op(sm, d[l]);
                        l++;
                    }
                }
                return l - size;
            }
            sm = op(sm, d[l]);
            l++;
        } while ((l & -l) != l);
        return _n;
    }
 
    template <bool (*g)(S)> int min_left(int r) {
        return min_left(r, [](S x) { return g(x); });
    }
    template <class G> int min_left(int r, G g) {
        assert(0 <= r && r <= _n);
        assert(g(e()));
        if (r == 0) return 0;
        r += size;
        for (int i = log; i >= 1; i--) push((r - 1) >> i);
        S sm = e();
        do {
            r--;
            while (r > 1 && (r % 2)) r >>= 1;
            if (!g(op(d[r], sm))) {
                while (r < size) {
                    push(r);
                    r = (2 * r + 1);
                    if (g(op(d[r], sm))) {
                        sm = op(d[r], sm);
                        r--;
                    }
                }
                return r + 1 - size;
            }
            sm = op(d[r], sm);
        } while ((r & -r) != r);
        return 0;
    }
 
  private:
    int _n, size, log;
    std::vector<S> d;
    std::vector<F> lz;
 
    void update(int k) { d[k] = op(d[2 * k], d[2 * k + 1]); }
    void all_apply(int k, F f) {
        d[k] = Map(f, d[k]);
        if (k < size) lz[k] = Com(f, lz[k]);
    }
    void push(int k) {
        all_apply(2 * k, lz[k]);
        all_apply(2 * k + 1, lz[k]);
        lz[k] = id();
    }
};
 
// [[ Assignment and Sum ]]   
using S = pll  ; using F = ll ;
S op( S a,S b) {return {a.ff+b.ff,a.ss+b.ss} ;}
S e() {return  {0,0}; }
S Map(F f, S x) { if(f==-1)return  x;
  else return {f*x.ss ,x.ss};}
F Com(F f, F  g) { return  f==-1?g:f ;}
F id() { return -1 ; }

// [[ Addition and Sum ]] 
using S = pll; using F = ll ;
S op( S a,S b) { return  {a.ff+b.ff,a.ss+b.ss};}
S e() { return  {0,0}; }
S Map(F f, S x) {return {x.ff + (x.ss*f ),x.ss};}
F Com(F f, F g) { return  f+g ; }
F id() { return 0; }


////


template <class S, S (*op)(S, S), S (*e)()> 
struct segtree {
    int n,size;
    vector< S > d;
    segtree(int n) : segtree(std::vector<S>(n, e())) {}
    segtree(const vector<S>& v){
        n=int(v.size()),size=1; 
        while(size<n)size<<=1;
        d=vector<S>(2 * size, e());
        for(int i=0;i<n;i++)d[size + i] = v[i];
        for(int i=size-1;i>=1;i--) update(i);
    }
    void update(int k) { 
        d[k] = op(d[2 * k], d[2 * k + 1]); 
    }
    void set(int p, S x) {
        assert(0 <= p && p < n);
        p += size;  d[p] = x;
        while(p > 1)update(p >>=1);
    }
    S get(int p) const {
        assert(0 <= p && p < n);
        return d[p + size];
    }

    S prod(int l, int r) const { // [l,r)
        assert(0 <= l && l <= r && r <= n);
         S sml = e(), smr = e();
          for(l += size ,r += size ; l<r ;  l>>=1,r>>=1 ){
              if (l & 1) sml = op(sml, d[l++]);
              if (r & 1) smr = op(d[--r], smr); 
         }
        return op(sml, smr);
    }
    S all_prod() const { return d[1]; }

    template <bool (*f)(S)> int max_right(int l) const {
        return max_right(l, [](S x) { return f(x); });
    }
    template <class F> int max_right(int l, F f) const {
        assert(0 <= l && l <= n);
        assert(f(e()));
        if (l == n) return n;
        l += size;
        S sm = e();
        do {
            while (l % 2 == 0) l >>= 1;
            if (!f(op(sm, d[l]))) {
                while (l < size) {
                    l = (2 * l);
                    if (f(op(sm, d[l]))) {
                        sm = op(sm, d[l]);
                        l++;
                    }
                }
                return l - size;
            }
            sm = op(sm, d[l]);
            l++;
        } while ((l & -l) != l);
        return  n;
    }

    template <bool (*f)(S)> int min_left(int r) const {
        return min_left(r, [](S x) { return f(x); });
    }
    template <class F> int min_left(int r, F f) const {
        assert(0 <= r && r <= n);
        assert(f(e()));
        if (r == 0) return 0;
        r += size;
        S sm = e();
        do {
            r--;
            while (r > 1 && (r % 2)) r >>= 1;
            if (!f(op(d[r], sm))) {
                while (r < size) {
                    r = (2 * r + 1);
                    if (f(op(d[r], sm))) {
                        sm = op(d[r], sm);
                        r--;
                    }
                }
                return r + 1 - size;
            }
            sm = op(d[r], sm);
        } while ((r & -r) != r);
        return 0;
    }


};

