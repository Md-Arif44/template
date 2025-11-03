 
template<class T> struct RMQ {
  V<V<T>> jmp;
  RMQ(const V<T>& V) : jmp(1, V) {
    for (int pw = 1, k = 1; pw * 2 <= sz(V); pw *= 2, ++k) {
      jmp.emplace_back(sz(V) - pw * 2 + 1);
      rep(j,sz(jmp[k]))
        jmp[k][j] = min(jmp[k - 1][j], jmp[k - 1][j + pw]);
    }
  }
  T query(int a, int b) { // [a,b)
    assert(a < b); // or return inf if a == b
    int dep = 31 - __builtin_clz(b - a);
    return min(jmp[dep][a], jmp[dep][b - (1 << dep)]);
} };
struct LCA {
  int T = 0;
  vi time, path, ret ,tour ;
  RMQ<int> rmq;
  LCA(vector<vi>& C) : time(sz(C)), rmq((dfs(C,0,-1), ret)) {}
  void dfs(vector<vi>& C, int v, int par) {
    time[v] = T++; tour.pb(v) ;
    for (int y : C[v]) if (y != par) {
      path.pb(v), ret.pb(time[v]);
      dfs(C, y, v);
    }
    tour.pb(v) ;
  }
 
  int lca(int a, int b) {
    if (a == b) return a;
    tie(a, b) = minmax(time[a], time[b]);
    return path[rmq.query(a, b)];
  }
  //dist(a,b){return depth[a] + depth[b] - 2*depth[lca(a,b)];}
};



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
