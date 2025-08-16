
// lde
template<typename T>
T extgcd(T a, T b, T &x, T &y)
{
    if (a == 0) return x = 0, y = 1, b;
    T p = b / a;
    T g = extgcd(b - p * a, a, y, x);
    x -= p * y;
    return g;
}
/// Return true if there exist such (x, y) satisfy ax + by = c
/// Find (&g) = gcd(a, b)
/// Find (&x, &y) satisfy ax + by = c
template<typename T>
bool find_any_solution(T a, T b, T c, T &x, T &y, T &g)
{
    if (a == 0 && b == 0) /// 0x + 0y = c
    {
        if (c != 0) return false;
        x = y = g = 0;
        return true;
    }

    if (a == 0) /// 0x + by = c
    {
        if (c % b != 0) return false;
        x = 0, y = c / b, g = abs(b);
        return true;
    }

    if (b == 0) /// ax + 0y = c
    {
        if (c % a != 0) return false;
        x = c / a, y = 0, g = abs(a);
        return true;
    }

    /// ax + by = c
    g = extgcd(abs(a), abs(b), x, y);
    if (c % g != 0) return false;

    x *= (a < 0 ? -1 : +1) * c / g;
    y *= (b < 0 ? -1 : +1) * c / g;
    return true;
}
/// Find the next/prev (cnt)-th solution of ax + by = c
template<typename T>
void shift_solution(T & x, T & y, T a, T b, T cnt)
{
    x += cnt * b;
    y -= cnt * a;
}
template<typename T = long long>
T find_all_solutions(T a, T b, T c, T min_x, T max_x, T min_y, T max_y) {
    if (min_x > max_x) return 0; /// Invalid range
    if (min_y > max_y) return 0; /// Invalid range

    if (a == 0 && b == 0) /// 0x + 0y = c
    {
        if (c != 0) return 0; /// No solution
        return 1LL * (max_x - min_x + 1) * (max_y - min_y + 1); /// Ways to select (x) and (y) in range
    }

    if (a == 0) /// 0x + by = c <=> y = c / b
    {
        if (c % b != 0) return 0; /// No solution
        if (1LL * min_y * b > c) return 0; /// Out of range: min > y
        if (1LL * max_y * b < c) return 0; /// Out of range: max < y
        return max_x - min_x + 1; /// Ways to select (x) in range    
    }

    if (b == 0) /// ax + 0y = c <=> x = c / a
    {
        if (c % a != 0) return 0; /// No solution
        if (1LL * min_x * a > c) return 0; /// Out of range: min > x
        if (1LL * max_x * a < c) return 0; /// Out of range: max < x
        return max_y - min_y + 1; /// Ways to select (y) in range    
    }

    T x, y, g;
    if (!find_any_solution(a, b, c, x, y, g)) return 0;
    a /= g;     
    b /= g;

    T sign_a = a > 0 ? +1 : -1;
    T sign_b = b > 0 ? +1 : -1;

    shift_solution(x, y, a, b, (min_x - x) / b);
    if (x < min_x) shift_solution(x, y, a, b, sign_b);
    if (x > max_x) return 0;
    T lx1 = x;

    shift_solution(x, y, a, b, (max_x - x) / b);
    if (x > max_x) shift_solution(x, y, a, b, -sign_b);
    T rx1 = x;

    shift_solution(x, y, a, b, -(min_y - y) / a);
    if (y < min_y) shift_solution(x, y, a, b, -sign_a);
    if (y > max_y) return 0;
    T lx2 = x;

    shift_solution(x, y, a, b, -(max_y - y) / a);
    if (y > max_y) shift_solution(x, y, a, b, sign_a);
    T rx2 = x;

    if (lx2 > rx2)swap(lx2, rx2);
    T lx = max(lx1, lx2);
    T rx = min(rx1, rx2);

    if (lx > rx) return 0;
    return (rx - lx) / abs(b) + 1;
}
Find the solution with minimum value x+y 
 x^=x+k (b/g);
  y^=y-k (a/g);
  if a<b select smalest possible k
  if b>a	select largest possible k
