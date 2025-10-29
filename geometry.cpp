template <class T> int sgn(T x) { return (x > 0) - (x < 0); }
const double eps = 1e-9;
const double PI = acos((double)-1.0);
// int sgn(double x) { return (x > eps) - (x < -eps); }

template<class T> struct Point { 
typedef Point P ;T x, y;
explicit Point(T _x=0, T _y=0) : x(_x),y(_y){}
bool operator<(P p) const { return tie(x,y) < tie(p.x,p.y); }
bool operator==(P p) const { return tie(x,y)==tie(p.x,p.y); }
P operator+(P p) const{return P(x+p.x,y+p.y);}
P operator-(P p) const{return P(x-p.x,y-p.y);}
P operator*(T d) const { return P(x*d, y*d); }
P operator/(T d) const { return P(x/d, y/d); }
T dot(P p) const { return x*p.x + y*p.y; }
T cross(P p) const { return x*p.y - y*p.x; }
T cross(P a, P b) const { return (a-*this).cross(b-*this); }
T dist2() const { return x*x + y*y; }
double dist() const { return sqrt((double)dist2()); }
// angle to x-axis in interval [-pi, pi]
double angle() const { return atan2(y, x); }
// makes dist() = 1
P unit() const { return *this/dist(); } 
// rotate by +90 degree
P perp() const { return P(-y, x); }
P normal() const { return perp().unit(); }
//rotate 'a' radians ccw around (0,0)
P rotate(double a) const { return P(x*cos(a)-y*sin(a),x*sin(a)+y*cos(a)); }
friend ostream& operator<<(ostream& os, P p) { 
  return os<<"("<< p.x << "," << p.y << ")";}
friend istream &operator>>(istream &in,P &p){return in >>p.x>>p.y; }
};
// a--b seg [cl,0] [-1, cw ] [1 ,ccw]
template<class P> 
int orientation(P a, P b, P c) { return sgn(a.cross(b,c)); }

double rad_to_deg(double r) { return (r * 180.0 / PI); }
double deg_to_rad(double d) { return (d * PI / 180.0); }

template<class P>  double get_angle(P a, P b) {
    double costheta = a.dot(b) / a.dist() / b.dist();
    return acos(max((double)-1.0, min((double)1.0, costheta)));
}
template<class P> 
bool is_point_in_angle(P b, P  a, P c, P  p) { // does point p lie in angle <bac
    assert(orientation(a, b, c) != 0);
    if (orientation(a, c, b) < 0) swap(b, c);
    return orientation(a, c, p) >= 0 && orientation(a, b, p) <= 0;
}
// checks is point p in s-- e line segment 
template<class P> 
bool onSegment(P s, P e, P p) {return p.cross(s, e) == 0 && (s - p).dot(e - p) <= 0;}
template<class P> 
vector<P> segInter(P a, P b, P c, P d) {
     auto oa = c.cross(d, a), ob = c.cross(d, b),
          oc = a.cross(b, c), od = a.cross(b, d);
     // Checks if intersection is single non-endpoint point.
     if (sgn(oa) * sgn(ob) < 0 && sgn(oc) * sgn(od) < 0)
          return {(a * ob - b * oa) / (ob - oa)};
     set<P> s;
     if (onSegment(c, d, a)) s.insert(a);
     if (onSegment(c, d, b)) s.insert(b);
     if (onSegment(a, b, c)) s.insert(c);
     if (onSegment(a, b, d)) s.insert(d);
     return {all(s)};
}
//Returns the shortest distance between point p and the line segment from point s to e.
// typedef Point<double> P;
// double segDist(P& s, P& e, P& p) {
//   if (s==e) return (p-s).dist();
//   auto d = (e-s).dist2(), t = min(d,max(.0,(p-s).dot(e-s)));
//   return ((p-s)*d-(e-s)*t).dist()/d;
// }


// [-1 for boundary ][1 for inside ] else outside
template<class P> 
int P_i_polygon(vector<P>&p ,P &x){
    P y = P(x.x + 1, ll(1e9)+7ll); 
    int cnt =0 ,n=sz(p) ;
    rep(i,n){
            cnt+=sz(segInter(p[i] ,p[(i+1)%n] ,x,y) )>0 ;
            if(onSegment(p[i], p[ (i+1)%n], x) )return -1 ;
     }
      if(cnt&1) return 1;
      else return 0 ;
}


template< class P>
vector<P> convex_hull(vector<P> &p) {
  if (p.size() <= 1) return p;
  vector<P> v = p;
    sort(v.begin(), v.end());
    vector<P> up, dn;
    for (auto& p : v) {
        while (up.size() > 1 && orientation(up[up.size() - 2], up.back(), p) >= 0) {
            up.pop_back();
        }
        while (dn.size() > 1 && orientation(dn[dn.size() - 2], dn.back(), p) <= 0) {
            dn.pop_back();
        }
        up.pb(p); dn.pb(p);
    }
    v = dn;
    if (v.size() > 1) v.pop_back();
    reverse(up.begin(), up.end());
    up.pop_back();
    for (auto& p : up)  v.pb(p);
    
    if (v.size() == 2 && v[0] == v[1]) v.pop_back();
    return v;
}
