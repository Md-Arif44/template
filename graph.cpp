
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
