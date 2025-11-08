/**
 * Author: Unknown
 * Date: 2002-09-15
 * Source: predates tinyKACTL
 * Description: Finds two integers $x$ and $y$, such that $ax+by=\gcd(a,b)$. If
 * you just need gcd, use the built in \texttt{\_\_gcd} instead.
 * If $a$ and $b$ are coprime, then $x$ is the inverse of $a \pmod{b}$.
 */

ll euclid(ll a, ll b, ll &x, ll &y) {
	if (!b) return x = 1, y = 0, a;
	ll d = euclid(b, a % b, y, x);
	return y -= a/b * x, d;
}

/**
 * Author: Simon Lindholm
 * Date: 2019-05-22
 * License: CC0
 * Description: Chinese Remainder Theorem.
 *
 * \texttt{crt(a, m, b, n)} computes $x$ such that $x\equiv a \pmod m$, $x\equiv b \pmod n$.
 * If $|a| < m$ and $|b| < n$, $x$ will obey $0 \le x < \text{lcm}(m, n)$.
 * Assumes $mn < 2^{62}$.
 * Time: $\log(n)$
 * Status: Works
 */
ll crt(ll a, ll m, ll b, ll n) {
	if (n > m) swap(a, b), swap(m, n);
	ll x, y, g = euclid(m, n, x, y);
	assert((a - b) % g == 0); // else no solution
	x = (b - a) % n * x % n / g * m + a;
	return x < 0 ? x + m*n/g : x;
}

/**
 * Author: Lukas Polacek
 * Date: 2009-09-28
 * License: CC0
 * Source: folklore
 * Description: Operators for modular arithmetic. You need to set {\tt mod} to
 * some number first and then you can use the structure.
 */


const ll mod = 1e9+7; // change to something else
struct Mint {
  ll x;
  Mint(ll xx) : x(xx) {}
  Mint operator+(Mint b) { return Mint((x + b.x) % mod); }
  Mint operator-(Mint b) { return Mint((x - b.x + mod) % mod); }
  Mint operator*(Mint b) { return Mint((x * b.x) % mod); }
  Mint operator/(Mint b) { return *this * invert(b); }
  Mint invert(Mint a) {
    ll x, y, g = euclid(a.x, mod, x, y);
    assert(g == 1); return Mint((x + mod) % mod);
  }
  Mint power(ll e) {
    if (!e) return Mint(1);
    Mint r = this-> power (e / 2); r = r * r;
    return e&1 ? *this * r : r;
  }
};
