

// lca
struct LCA{
  int n ,Log;
  vector<int> dis ;
  vector<vector<int>>g,up;
  LCA(int n):n(n){
      Log= 32 - __builtin_clz(n); 
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
// for dfs order scc etc 
void  DFS(const vvi &G ,const vi&order){
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
    
  return ;
}
