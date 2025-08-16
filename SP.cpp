
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
  int n,m;
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


