
ll euclid(ll a, ll b, ll &x, ll &y) {
  if (!b) return x = 1, y = 0, a;
  ll d = euclid(b, a % b, y, x);
  return y -= a/b * x, d;
}
/// Return true if there exist such (x, y) satisfy ax + by = c
/// Find (&g) = gcd(a, b)
/// Find (&x, &y) satisfy ax + by = c
// for many solution  (x0 = x0+ k*b/g ) and (y0 = y0- k* a/g)  k is  a integer 
bool find_any_solution(int a, int b, int c, int &x0, int &y0, int &g) {
    g = euclid (abs(a), abs(b), x0, y0);
    if (c % g) return false;
    
    x0 *= c / g;
    y0 *= c / g;
    if (a < 0) x0 = -x0;
    if (b < 0) y0 = -y0;
    return true;
}
