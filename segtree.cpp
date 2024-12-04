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
