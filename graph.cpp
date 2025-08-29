
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
      if (depth[v] == -1) {
          ++attempt;
          depth[v] = 0;
          root[v] = v;
          pa[v] = -1;
          dfs(v);
      
      }
    }
    assert((int) order.size() == n);
  }
  void add(int u,int v){
         assert(u<n && v<n && u>=0 && v>=0);
         g[u].push_back(v);
  }

};


// dfs 
struct Graph{
      int n,attempt=0;
      vector<vector<int>> g;
      vector<int>pa,was,pos,end,order,root,sz,dist;
     Graph (int _n):n(_n) {
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
  }
  void dfs_v(int v){
       if (pos[v] == -1) {
          ++attempt;
          dist[v] = 0;
          root[v] = v;
          pa[v] = -1;
          dfs(v);
      }
  }
  void dfs_all() {
    for (int v = 0; v < n; v++) {
      if (pos[v] == -1) {
          ++attempt;
          dist[v] = 0;
          root[v] = v;
          pa[v] = -1;
          dfs(v);
      
      }
    }
    assert((int) order.size() == n);
  }
  void add(int u,int v){
         assert(u<n && v<n && u>=0 && v>=0);
         g[u].push_back(v);
  }
 
};
