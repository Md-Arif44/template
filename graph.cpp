
art bri     
struct graph{
      int n,attempt;
      vector<vector<int>> g;
      vector<int>pa,min_depth,depth; 
      vector<int>was,order,root,sz,dist;
     graph (int _n):n(_n) {
           assert(n>0); 
           g=vector<vector<int>>(n);        
           pa= vector<int>(n, -1);
           order.clear();
           sz = vector<int>(n, 0);
           root = vector<int>(n, -1);
           depth = vector<int>(n, -1);
           min_depth = vector<int>(n, -1);
           dist = vector<int>(n);
           was = vector<int>(n, -1);
           attempt = 0;
     }

  void is_bridge(int fo,int to){
  }
  void is_artpoint(int fo){
  }
  void dfs(int v){
    was[v] = attempt; 
    order.push_back(v);
    sz[v] = 1;
    min_depth[v] = depth[v];
    int child=0;
    for (auto to : g[v]) {
    
      if (to == pa[v]) {
          continue;
      }
      if (was[to] == attempt) {
           min_depth[v] = min(min_depth[v], depth[to]);       
          continue;
      }
      depth[to] = depth[v] + 1;
      pa[to] = v;
      root[to]=(root[v] != -1 ? root[v] : to);
      dfs(to);
      sz[v] += sz[to];
      min_depth[v] = min(min_depth[v], min_depth[to]);
      
        if(min_depth[to]>=depth[v]){
           is_artpoint(v);
        }
        if(min_depth[to]>depth[v]){
           is_bridge(v,to);
        }   
        ++child;
    }
    
     
      if(child>1 && pa[v]==-1){ 
         is_artpoint(v);
      }   
  }
  void dfs_v(int v){
       if (depth[v] == -1) {
         ++attempt;
         depth[v] = 0;
         root[v] = v;
          pa[v] = -1;
          dfs(v);
      }
  }
  void dfs_all() {
    for (int v = 0; v < n; v++) {
       dfs_v(v);
    }
    assert((int) order.size() == n);
  }
  void add(int u,int v){
         assert(u<n && v<n && u>=0 && v>=0);
         g[u].push_back(v);
  }

};
// dfs scc 
 
struct Graph{
      int n,attempt=0;
      vector<vector<int>> g;
      vector<int>pa,was,pos,end,order,end_order, root,sz,dist;
     Graph(int _n):n(_n) {
           assert(n>0); 
           g.assign(n,{});      
           pa.resize(n, -1);
           pos.resize(n,-1);
           end.resize(n,-1);
           sz.resize(n, 0);
           root.resize(n, -1);
           dist.resize(n);
           was.resize(n, -1);
     }
  void dfs(int v){
    was[v] = attempt; 
     pos[v]=int(order.size()) ;
         order.push_back(v);
         sz[v]=1;
         for(auto to: g[v]){
             if(pos[to]==-1){
                pa[to]=v;
                dist[to]=dist[v]+1;
                root[to]=(root[v]!= -1?root[v] : to);
                dfs(to);
                sz[v]+=sz[to];
             }
         }
      end[v]=int(order.size())-1 ;
      end_order.pb(v);
  }
  void dfs_v(int v){
       if (pos[v] == -1) {
          ++attempt;dist[v] = 0;
          root[v] = v;pa[v] = -1;
          dfs(v);
      }
  }
  void dfs_all() {
    for (int v = 0; v < n; v++) {
         dfs_v(v);
    }
    assert((int) order.size() == n);
  }
  void add(int u,int v){
      assert(u<n && v<n && u>=0 && v>=0);
         g[u].push_back(v);
  }
 
};
 
 struct Graph{
      int n,attempt=0;
      vector<vector<int>> g;
      vector<int>pa,was,pos,end,order,end_order, root,sz,dist;
     Graph(int _n):n(_n) {
           Resize(); 
           g.assign(n,{});      
      }
     void Resize(){
           assert(n>0);
           pa.resize(n, -1);
           pos.resize(n,-1);
           end.resize(n,-1);
           sz.resize(n, 0);
           root.resize(n, -1);
           dist.resize(n);
           was.resize(n, -1);
    
    }
    void Clear(){
          pa.clear();attempt=0;
           pos.clear();
           end.clear();
           sz.clear();
           root.clear();
           dist.clear();
           was.clear();
          order.clear();
          end_order.clear();
    }
    void dfs(int v){
     was[v] = attempt; 
      pos[v]=int(order.size()) ;
         order.push_back(v);
         sz[v]=1;
         for(auto to: g[v]){
             if(pos[to]==-1){
                pa[to]=v;
                dist[to]=dist[v]+1;
                root[to]=(root[v]!= -1?root[v] : to);
                dfs(to);
                sz[v]+=sz[to];
             }
         }
      end[v]=int(order.size())-1 ;
      end_order.pb(v);
  }
  void dfs_v(int v){
       if (pos[v] == -1) {
          ++attempt;dist[v] = 0;
          root[v] = v;pa[v] = -1;
          dfs(v);
      }
  }
  void dfs_all() {
    for (int v = 0; v < n; v++) {
         dfs_v(v);
    }
    assert((int) order.size() == n);
  }
  void add(int u,int v){
      assert(u<n && v<n && u>=0 && v>=0);
         g[u].push_back(v);
  }
 
};
 




// lca
struct LCA{
  int n ,Log;
  vector<int> dis ;
  vector<vector<int>>g,up;
  LCA(int n):n(n){
      Log=Msb(n) +1;
      up.resize(n,vector<int>(Log)); 
      g.assign(n,{}); dis.resize(n);
  }
  void dfs(int v,int pa){
           up[v][0]=pa; 
           for(int i=1;i<Log;i++){
               if(dis[v]>= (1<<i)) up[v][i]=up[up[v][i-1]][i-1];
           }
           for(auto u: g[v]){
             if(u==pa)continue ;
               dis[u]=dis[v]+1;
               dfs(u,v) ;
           }

  }
  int lca(int u,int v){
           if(dis[u]<dis[v])swap(u,v);
 
           int k=dis[u]-dis[v] ;
          for(int i=0;i<Log;i++){
             if(k&(1<<i)){
                 u=up[u][i];
             }
          }
         if(u==v)return u ;
          for(int i=Log-1;i>=0 ;i-- ){
             if(up[u][i]!=up[v][i]){
                u=up[u][i];
                v=up[v][i];
             }
          }
         return up[v][0];
 
  }
  void build(int root ){
    dfs(root,0);
  }
  void add(int u,int v){
      assert(u<n && v<n && u>=0 && v>=0);
         g[u].push_back(v);
  }

};

pair<vi,vi>  DFS(const vvi &G ,const vi&order){
   int n=sz(G),attempt=0 ;
   vi pa(n,-1),was(n,-1),pos(n,-1),end(n,-1);
   vi ord,root(n,-1),Sz(n),dist(n),e_ord;   

   auto dfs=[&] (auto &&dfs,int v)-> void {
      was[v] = attempt; 
       pos[v]=sz(ord) ;
         ord.pb(v);
         Sz[v]=1;
         for(auto to: G[v]){
             if(pos[to]==-1){
                pa[to]=v;
                dist[to]=dist[v]+1;
                root[to]=(root[v]!= -1?root[v] : to);
                dfs(dfs, to);
                Sz[v]+=Sz[to];
             }
         }
      end[v]=sz(ord)-1 ;
      e_ord.pb(v);
  };
     rep(i,n){
        int v=(order[i]==-1? i:order[i]) ;
         if (pos[v] == -1) {
          ++attempt;dist[v] = 0;
          root[v] = v;pa[v] = -1;
          dfs(dfs, v);
        }
     }
    reverse(all(e_ord));
  return make_pair(e_ord,was);
}
