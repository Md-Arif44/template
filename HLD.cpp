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
