 bfs 
 vector<int>pa(n,-1),dis(n,-1),vis(n);
 auto bfs=[&] (int s)-> void {
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
};
 2D
 int n,m;
       vector<vector<string>>  girid(n);
       vector level (n, vector<int>(m));
       map<pair<int,int> ,pair<int,int>>parant;
       vector<char>Direction;
      
       auto is_valid=[&](int x,int y)-> bool { return x<n && y<m && girid[x][y]!='#'   ;};
   
        auto bfs=[&](pair<int ,int> s)-> void {


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
                    if(is_valid(x,y) && level[x][y]==-1){
                        q.push({x,y});
                        level[x][y]=level[u.ff][u.ss]+1;
                        parant[{x,y}]={u.ff,u.ss};
                    }
                }
           }
      };

Dijkstra

     struct edge{
         ll fo,w;
         edge(ll fo,ll w):fo(fo),w(w){}
         bool operator<(edge const &e)const {
          return e.w<w;
          }
        };
        int n,m;
        cin>>n>>m;
        vector<vector< edge >>g(n);
        vector<ll>dis(n,inf);
        int u,v;
        ll w;
        for(int i=0;i<m;i++){ 
          cin>>u>>v>>w;
          --u,--v;
         g[u].push_back({v,w});
       
        }
        priority_queue<edge>q;
        q.emplace(0,0);
        dis[0]=0;
         while(!q.empty()){
           auto v=q.top();
           q.pop();
           if(v.w>dis[v.fo])continue ;    
          
           for(auto &[u,w]:g[v.fo]){
                if(w+dis[v.fo]<dis[u]){ 
                  dis[u]=w+dis[v.fo];
                  q.emplace(u,dis[u]);
                }
           }          
        }
     
 Floyd Warshall
  int n,m,q;
  vector<vector<ll>>dis(n,vector<ll>(n,inf));
  
  for(int i=0;i<n;i++)dis[i][i]=0;
    
  for(int k=0;k<n;k++){
     for(int i=0;i<n;i++){
       for(int j=0;j<n;j++){
          if(dis[i][k]!=inf && dis[k][j]!=inf ){
            dis[i][j]=min(dis[i][j],dis[i][k]+dis[k][j]);
          }
       }
     }
  }


