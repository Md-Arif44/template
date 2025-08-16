
template <class S, S (*op)(S, S), S (*e)()> 
struct segtree {
  public:
    segtree() : segtree(0) {}
    explicit segtree(int n) : segtree(std::vector<S>(n, e())) {}
    explicit segtree(const std::vector<S>& v) : _n(int(v.size())) {
      if(_n&(_n-1)== 0)size=1<<Msb(_n);
      else size=1<<(Msb(_n)+1);
      
        log =__builtin_ctz((unsigned int)size);
        d = std::vector<S>(2 * size, e());
        for (int i = 0; i < _n; i++) d[size + i] = v[i];
        for (int i = size - 1; i >= 1; i--) {
            update(i);
        }
    }

    void set(int p, S x) {
        assert(0 <= p && p < _n);
        p += size;
        d[p] = x;
        for (int i = 1; i <= log; i++) update(p >> i);
    }

    S get(int p) const {
        assert(0 <= p && p < _n);
        return d[p + size];
    }

    S prod(int l, int r) const {
        assert(0 <= l && l <= r && r <= _n);
        S sml = e(), smr = e();
        l += size;
        r += size;

        while (l < r) {
            if (l & 1) sml = op(sml, d[l++]);
            if (r & 1) smr = op(d[--r], smr);
            l >>= 1;
            r >>= 1;
        }
        return op(sml, smr);
    }

    S all_prod() const { return d[1]; }

    template <bool (*f)(S)> int max_right(int l) const {
        return max_right(l, [](S x) { return f(x); });
    }
    template <class F> int max_right(int l, F f) const {
        assert(0 <= l && l <= _n);
        assert(f(e()));
        if (l == _n) return _n;
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
        return _n;
    }

    template <bool (*f)(S)> int min_left(int r) const {
        return min_left(r, [](S x) { return f(x); });
    }
    template <class F> int min_left(int r, F f) const {
        assert(0 <= r && r <= _n);
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

  private:
    int _n, size, log;
    std::vector<S> d;
    void update(int k) { d[k] = op(d[2 * k], d[2 * k + 1]); }
};
// l 0 base r 1 base 
int op(int x ,int y){
   return x+y;
}
int e(){
   return 0 ;
}





struct Node{
     
};
struct segtree{
  vector<Node>tree;
   int size;
   Node NUTURL_ELEMENT={};
   ll no_oparation ;
     segtree(int n){
           assert(n>0);
            size=1;
           while(size<n)size<<=1;
            tree.assign(size<<1,{});
            size=n-1;
        }
        template <class T>     
        inline void leaf(Node &a,T value){

        }
        template<class T>
        inline void apply(Node &a,T value,int lx=0,int rx=0){
           
        }
        inline void push(int x,int lx,int rx){

        /*
           if(lx==rx)return;
           if(tree[x].add!=no_oparation){
             apply(tree[(x<<1)+1],tree[x].add);
             apply(tree[(x<<1)+2],tree[x].add);
             tree[x].add=no_oparation;
           }

          */ 
        }
       inline Node unite(const Node &a,const  Node &b){
         
          
        }
        template<class T>
       void build (vector<T>&v,int x, int lx,int rx){
            if(lx==rx && lx<sz(v)){
           leaf(tree[x],v[lx]);
            return;
            }
        int mid =(lx+rx)>>1;
            build(v,(x<<1)|1,lx ,mid);
            build (v,(x<<1)+2,mid+1,rx);

        tree[x]=unite(tree[(x<<1)|1],tree[(x<<1)+2]);
        }
        template<class T>
        void build(vector<T>&a){
            build(a,0,0,size);
        }
        template<class T>
        void set(int x , int lx ,int rx, int l,int r, T &value){
           push(x,lx,rx);
            if (l > rx || r < lx)
                return ; 
            if (lx >= l && rx <= r){
               apply(tree[x],value,lx,rx);         
                return;
            }  
            int mid =(lx+rx)>>1;
           set((x<<1)|1,lx,mid,l,r,value);
            set((x<<1)+2,mid+1,rx,l,r,value);
            tree[x]=unite(tree[(x<<1)|1],tree[(x<<1)+2]);
        }
        template<class T>
        void set(int l ,int r,T &val){
          assert(0 <= l && l <= r && r <= size);
            set(0,0,size,l,r,val);
        }
        Node get(int x, int lx, int rx, int l, int r){
            push(x,lx,rx);
            if (l > rx || r < lx)
                return NUTURL_ELEMENT; 
            if (lx >= l && rx <= r)
                return tree[x];

            int mid =(lx+rx)>>1;
          Node p1 = get((x<<1)|1,lx,mid,l,r);
          Node  p2 = get((x<<1)+2,mid+1,rx,l,r);
            return unite(p1,p2); 
        }
       Node get(int l,int r){
        assert(0 <= l && l <= r && r <= size);
           return get(0,0,size,l,r);
        }
        
};


