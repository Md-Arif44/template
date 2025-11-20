
void MO_qary(){
   
/*
  mo tree quary  make 2n tour vec 
  for path u to v if  u == lca(u,v) then quary = in[u] to in [v ] ;
  else quary  out[u] to in[v] then add lca node  
  all are 0 base 
*/

const int Block = 335  ;
const int N =  100009;

struct  MO{
  int l,r , id ;
  bool operator<(const MO &e )const {
   return  make_pair(l/Block , r)< make_pair(e.l/Block ,e.r) ;
  }
} quary[N];
/*
void add (int id){   }
void rem(int id){ }
int clc(){ }
*/
vi tour , in(n),out(n) ;
 
    auto dfs=[&](auto dfs,int v,int pa)-> void {
        in[v]=sz(tour);
        tour.pb(v) ;
        for(auto u: G[v]){
           if(u!= pa ){
             dfs(dfs,u,v) ;
           }
        }
       out[v]=sz(tour);
       tour.pb(v);
    }; dfs(dfs,0,-1);
    
  

 rep(i,q ){ // number of quary 
  int u ,v ;
   cin>>u>>v ;
   --u,--v;
    if(in[u]> in[v] )swap(u,v) ;
      int x= A.lca(u,v) ;
    if( x == u){
       quary[i]={in[u],in[v],i};
    }else{
       quary[i]={out[u],in[v],i};
    }
 }
 sort(quary,quary+q);
  
int cl  = 0 ,cr=  -1 ;
vi ans (q) ;
for (int i=0;i<q;i++) {
  auto &[l,r ,k, lca, id ]  = quary[i];     
    while (cl > l)add( --cl);
    while (cr < r)add(++cr );
    while (cl < l)rem(cl++ );
    while (cr > r)rem (cr-- );
    ans[id]= clc () ;
}


}
 

void MO_qary(){
  // mo  quary 0 base
const int Block = 335 ,N =  100009;

struct  MO{
   int l,r , id ;
   bool operator<(const MO &e )const {
     return  make_pair(l/Block , r)< make_pair(e.l/Block ,e.r) ;
    }
} quary[N];
/*
void add (int id){   }
void rem(int id){ }
int clc(){ }
*/
rep(i,q ){ // number of quary 
 
   int l,r ;
   cin>>l>>r;
   --l,--r; 
   quary[i]={l,r,i};
 }
 sort(quary,quary+q);
 int cl  = 0 ,cr=  -1 ;
vi ans (q) ;
for (int i=0;i<q;i++) {
  auto &[l,r ,k, lca, id ]  = quary[i];     
    while (cl > l)add( --cl);
    while (cr < r)add(++cr );
    while (cl < l)rem(cl++ );
    while (cr > r)rem (cr-- );
    ans[id]= clc () ;
}


}
