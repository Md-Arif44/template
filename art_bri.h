
int n , timer;
V<vi> adj; 
V<pii > bri ;
vi vis, tin,Sz, low ,art ;
void dfs(int v, int p = -1) {
    vis[v] = true;
    Sz[v]=1  ;
    tin[v] = low[v] = timer++;
    int child =0 ;
    bool pas = false;
    for (int to : adj[v]) {
        if (to == p && !pas) {
            pas = true; 
            continue;
        }
        if (vis[to]) {
            low[v] = min(low[v], tin[to]);
        } else {
            dfs(to, v);
            Sz[v]+=Sz[to] ;
            low[v] = min(low[v], low[to]);
            if (low[to] >= tin[v] && p!=-1)art[v]=1;
            if (low[to] > tin[v])bri.pb({v, to});
             ++child;
        }
    }
    if(p == -1 && child > 1)art[v ] =1 ;
}

void art_bri() {

    timer = 0; bri.clear() ;
    Sz=vi(n,0) ;
    vis=vi(n, 0);art=vi(n,0);
    tin= vi (n, -1);low=vi(n, -1);
    rep (i,n) {
      if(!vis[i])dfs(i);
    }

}
