// ord =  ord[pos[i]] = a[i], ord for segtree   ; 
struct HLD {
  V<V<int>> adj ;
  vi pa, dth, havy, head, pos;
  int cur_pos,n;

int dfs(int v ) {
    int size = 1;
    int max_c_size = 0;
    for (int c : adj[v]) {
        if (c != pa[v]) {
            pa[c] = v, dth[c] = dth[v] + 1;
            int c_size = dfs(c);
            size += c_size;
            if (c_size > max_c_size)
                max_c_size = c_size, havy[v] = c;
        }
    }
 return size;
}

void decompose(int v, int h) {
    head[v] = h, pos[v] = cur_pos++;
    if (havy[v] != -1) decompose(havy[v], h);
    for (int c : adj[v]) {
        if (c != pa[v] && c != havy[v]) decompose(c, c);
    }
}

 HLD (V<V<int>> &G ): n(sz(G)){
    adj=G; cur_pos = 0;
    pa=vi(n),dth=vi(n),havy=vi(n,-1),head=vi(n),pos=vi(n);
    dfs(0);  decompose(0, 0);
}
// path [a,b]
 V<pii> path (int a, int b) {
  V< pii > Seg ;
    for (; head[a] != head[b]; b = pa[head[b]]) {

          if(dth[head[a]] > dth[head[b]]) swap(a, b);
      
      Seg.pb ({pos[head[b]], pos[b]});
    }
      if(dth[a] > dth[b]) swap(a, b);
    
      Seg.pb({pos[a], pos[b]});
    
    return Seg;
  }
} ;


// build  ord[ Id[i] ] =  a[i] for  tree 
template <bool VALS_EDGES = false > struct HLD {
  int N, tim = 0;
  vector<vi> adj;
  vi par, siz, rt, pos;

  HLD(vector<vi> adj_)
    : N(sz(adj_)), adj(adj_), par(N, -1), siz(N, 1),
      rt(N),pos(N) { dfsSz(0); dfsHld(0); }
  void dfsSz(int v) {
    for (int& u : adj[v]) {
      adj[u].erase(find(all(adj[u]), v));
      par[u] = v;
      dfsSz(u);
      siz[v] += siz[u];
      if (siz[u] > siz[adj[v][0]]) swap(u, adj[v][0]);
    }
  }
  void dfsHld(int v) {
    pos[v] = tim++;
    for (int u : adj[v]) {
      rt[u] = (u == adj[v][0] ? rt[v] : u);
      dfsHld(u);
    }
  }
  void process(int u, int v, vii &op) {
    for (;; v = par[rt[v]]) {
      if (pos[u] > pos[v]) swap(u, v);
      if (rt[u] == rt[v]) break;
      op.pb({pos[rt[v]], pos[v] + 1});
    }
     op.pb({pos[u] + VALS_EDGES, pos[v] + 1});
  }
  vii Path(int u, int v) { // return vec of [l,r)  
   vii op ; process(u, v,op) ;
    return op ;
  }
  int Id (int v){ //  return i such ord[i ] = a[i]  for tree 
     return pos[v] ;
  }
  pii querySubtree(int v) { // return [l,r)    
    return  {pos[v] + VALS_EDGES, pos[v] + siz[v]} ;
  }
};



