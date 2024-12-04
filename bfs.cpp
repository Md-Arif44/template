2D girid
/*
int n,m;
string girid[1005];
int level [1005][1005];
int dx[]={1,-1,0,0};
int dy []={0,0,1,-1};
map<pair<int,int> ,pair<int,int>>parant;
vector<char>Direction;
bool isvalid(int x, int y){
     return x>=0 and x<n and y>=0 and y<m and girid[x][y]!='#';
 }
void bfs(pair<int ,int> s){
    queue<pair<int,int>>q;
     q.push(s);
     level[s.ff][s.ss]=0;
     parant[{s.ff,s.ss}]={-1,-1};
     while(!q.empty()){
      pii u=q.front();
      q.pop();
          for(int i= 0;i<4;i++){
            int x=u.ff+dx[i];
            int y=u.ss+dy[i];
              if(isvalid(x,y) && level[x][y]==-1){
                  q.push({x,y});
                  level[x][y]=level[u.ff][u.ss]+1;
                  parant[{x,y}]={u.ff,u.ss};
              }
          }
     }
     return ;
}


*/



/*
struct graph{
      int n,m,g_time=0;
      vector<vector<int>>g;
      vector<bool>vis;
      vector<int>pa,dis,col,order; 
     void resize(){
          assert(n>0);
          g.assign(n+1,{});
          dis.resize(n+1);
          pa.resize(n+1,-1);
          order.clear();
          col.resize(n+1);
          vis.resize(n+1);
         assert(pa.size()>0 && dis.size()>0 && col.size()>0 && vis.size()>0);
     }
     void dfs(int s){
          order.emplace_back(s);
              vis[s]=1;
              for(auto i:g[s]){
                if(!vis[i]){
                   pa[i]=s;
                   dis[i]=dis[s]+1;
                   dfs(i);
                }
              }
     }
     void dfs_all(){
       for(int i=1;i<=n;i++){ 
         if(!vis[i]){
           dfs(i);
         }
       }
       assert(order.size()==n);
     }
     void bfs(int s){
            queue<int> q;
            q.emplace(s);
            vis[s] = true;
            pa[s] = -1;
             dis[s]=0 ;
            while (!q.empty()) {
                int v = q.front();
                q.pop();
                for (int u : g[v]) {
                    if (!vis[u]) {
                        vis[u] = true;
                        q.emplace(u);
                        dis[u] = dis[v] + 1;
                        pa[u] = v;
                    }
                }
            }
        }
  void Test_case(){
  
  }
};
*/
/*
  struct edge{
      int fo,w;
      edge(int fo,int w):fo(fo),w(w){}
      bool operator<(edge const &e)const {
        return e.w<w;
      }
};
void  Dijkstra(){
        int n,m;
         cin>>n>>m;
        vector<vector<pair<int,int>>>g(n+1);
        vector<int>dis(n+1,INT32_MAX),pa(n+1);
        
         for(int i=0;i<m;i++){ 
         int u,v,w;cin>>u>>v>>w;
         g[u].emplace_back(v,w);
         g[v].emplace_back(u,w);
        }
        priority_queue<edge>q;
        int sur=1;
        q.emplace(sur,0);
        dis[sur]=0;
         while(!q.empty()){
           auto v=q.top();
           q.pop();
           for(auto &[u,w]:g[v.fo]){
                if(w+dis[v.fo]<dis[u]){ 
                  dis[u]=w+dis[v.fo];
                  q.emplace(u,dis[u]);
                }
           }          
        }
}
*/
