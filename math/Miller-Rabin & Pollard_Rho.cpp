#include <bits/stdc++.h>

using namespace std;
typedef long long ll;

const int Times = 10;
const int N = 5500;

ll ct;
ll fac[N];

ll gcd(ll a, ll b) {
    return b ? gcd(b, a % b) : a;
}

ll multi(ll a, ll b, ll m) {
    ll ret = 0;
    a %= m;
    while (b) {
        if (b & 1) {
            ret = (ret + a) % m;
        }
        b >>= 1;
        a = (a + a) % m;
    }
    return ret;
}

ll qpow(ll a, ll b, ll m) {
    ll ret = 1;
    a %= m;
    while (b) {
        if (b & 1) {
            ret = ret * a % m;
        }
        b >>= 1;
        a = a * a % m;
    }
    return ret;
}

bool Miller_Rabin(ll n) {
    if (n == 2) return true;
    if (n < 2 || !(n & 1)) return false;
    ll m = n - 1;
    int k = 0;
    while ((m & 1) == 0) {
        k++;
        m >>= 1;
    }
    for (int i = 0; i < Times; ++i) {
        ll a = rand() % (n - 1) + 1;
        ll x = qpow(a, m, n);
        ll y = 0;
        for (int j = 0; j < k; ++j) {
            y = multi(x, x, n);
            if (y == 1 && x != 1 && x != n - 1) return false;
            x = y;
        }
        if (y != 1) return false;
    }
    return true;
}

ll pollard_rho(ll n, ll c) {
    ll i = 1, k = 2;
    ll x = rand() % (n - 1) + 1;
    ll y = x;
    while (true) {
        i++;
        x = (multi(x, x, n) + c) % n;
        ll d = gcd((y - x + n) % n, n);
        if (1 < d && d < n) return d;
        if (y == x) return n;
        if (i == k) {
            y = x;
            k <<= 1;
        }
    }
}

void find(ll n, int c) {
    if (n == 1) return;
    if (Miller_Rabin(n)) {
        fac[ct++] = n;
        return;
    }
    ll p = n;
    ll k = c;
    while (p >= n) p = pollard_rho(p, c--);
    find(p, k);
    find(n / p, k);
}

int main() {
    ll n;
    cin >> n;
    find(n, 120);
    sort(fac, fac + ct);
    //排好序的所有质因子 例如60被拆解为 2 2 3 5
}