struct graph{
      int n,attempt,g_time=0;
      vector<vector<int>>g;
      vector<int>pa,i_time,o_time,min_depth,depth,eular_depth; 
      vector<int>was,order,root,sz,pos,end,dist,eular;
     void resize(){
           assert(n>0); 
           g=vector<vector<int>>(n+1);        
           i_time=vector<int>(n+1,-1);
           o_time=vector<int>(n+1,-1);
           pa= vector<int>(n+1, -1);
           order.clear();
           pos = vector<int>(n+1, -1);
           end = vector<int>(n+1, -1);
           sz = vector<int>(n+1, 0);
           root = vector<int>(n+1, -1);
           depth = vector<int>(n+1, -1);
           min_depth = vector<int>(n+1, -1);
           dist = vector<int>(n+1);
           was = vector<int>(n+1, -1);
           attempt = 0;
     }
  void clear() {
        pa.clear();
        order.clear();
        pos.clear();
        end.clear();
        sz.clear();
        root.clear();
        depth.clear();
        min_depth.clear();
        dist.clear();
        was.clear();
        eular.clear();
        eular_depth.clear();
        i_time.clear();
        o_time.clear();
  } 

  void is_bridge(int fo,int to){
  }
  void is_artpoint(int fo){
  }
  void dfs(int v){
    was[v] = attempt; 
    pos[v] = (int) order.size();
    i_time[v]=g_time++;
    order.push_back(v);
    eular.push_back(v);
    eular_depth.push_back(depth[v]);
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
      eular.push_back(v);
      eular_depth.push_back(depth[v]);
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
     o_time[v]=g_time++;
     end[v] = (int) order.size() - 1;        
     if(child>1 && pa[v]==-1){ 
       is_artpoint(v);
     }   
  }
  void dfs_all() {
    for (int v = 1; v <= n; v++) {
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
  void Test_case(){ 
     
  }
};
