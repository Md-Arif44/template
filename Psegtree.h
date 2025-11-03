struct node{
  int64_t sum ;
  node *lc; node *rc;
  node(){ rc=lc=NULL;sum =0 ;
  }
};
node* build(int L,int R){
  node * cur  = new node() ;      
   if(L == R){
      cur->sum=a[L ];
     return  cur;
   }
   int m=(L+R)/2;
   cur->lc = build(L,m);
   cur->rc= build(m+1,R);
   cur->sum = cur->lc-> sum + cur->rc->sum ;  
   return cur ;  
}
ll quary (node *cur, int a,int b, int L,int R){
     
    if(a>R  || b< L )return 0ll ;
    if(L>= a && R<=b )return cur-> sum ;
   int m=(L+R)/2;
   ll p1 =quary(cur->lc ,a,b ,L,m);
   ll q1 =quary( cur->rc,a,b ,m+1,R);
  return p1+ q1 ;
}
node* update(node *cur ,int i,int val ,int L,int R){
   if(L==R ){
     node * ne = new node() ;
     ne ->sum= val ;
     return  ne ;
   }
   int m=(L+R)/2;
   if( i<= m ){
    node * ne  = new node() ;
    ne-> lc =update(cur->lc, i,val ,L,m) ;
    ne-> rc = cur-> rc ;
    ne-> sum = ne->lc-> sum + ne->rc->sum ;  
    return ne ;
  }else{
    node * ne  = new node() ;
    ne-> rc =update(cur->rc, i,val ,m+1,R ) ;
  
    ne-> lc = cur-> lc ;
    ne-> sum = ne->lc-> sum + ne->rc->sum ;  
     return ne ; 
   }
  
}
