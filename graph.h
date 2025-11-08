
pair<vi,vi> Fun (const V<vi> &g ,const vi &order  ){
     vi topo ,vis(sz(g)),was(sz(g))  ;
     int attempt=0;
     auto dfs=[&](auto dfs,int v)-> void {
             vis[v]=1;
             was[v]=attempt ;
            for(auto u: g[v]){
               if(vis[u]==0)
                  dfs(dfs, u);
            }
        topo.pb(v);
     };
     rep(i ,sz(order) ) {
           if(vis[order[i]] ==0 ){
                attempt++ ;
              dfs(dfs, order[i] ) ;
           }
     }
    reverse(all(topo))  ;
    return make_pair(topo,was) ;
}
// return SCC 
vi SCC(const V<vi> &g ,V <vi>&rg){
    vi ord(sz(g)) ; iota(all(ord),0) ;
    auto [topo,_ ]= Fun(g, ord);
    auto [__,was]= Fun(rg,topo); // was has scc with same val 
  return was; 
}


/**
 * Author: Unknown
 * Date: 2002-09-13
 * Source: predates tinyKACTL
 * Description: Topological sorting. Given is an oriented graph.
 * Output is an ordering of vertices, such that there are edges only from left to right.
 * If there are cycles, the returned list will have size smaller than $n$ -- nodes reachable
 * from cycles will not be returned.
 * Time: $O(|V|+|E|)$
 * Status: stress-tested
 */
 
vi topoSort(const vector<vi>& gr) {
  vi indeg(sz(gr)), q;
  for (auto& li : gr) for (int x : li) indeg[x]++;
  rep(i,sz(gr)) if (indeg[i] == 0) q.pb(i);
  rep(j,sz(q)) for (int x : gr[q[j]])
    if (--indeg[x] == 0) q.pb(x);
  return q;
}
