/* exBSGS算法
 * N ^ x = M (mod P), 其中N, P互质
 * 返回的x为最小解
 */
#include<bits/stdc++.h>

using namespace std;
typedef long long ll;
unordered_map<ll, int> H;

ll gcd(ll a, ll b) {
    if (!b) return a;
    return gcd(b, a % b);
}

ll expow(ll a, ll b, ll mod) {
    ll res = 1;
    while (b) res = ((b & 1) ? res * a % mod : res), a = a * a % mod, b >>= 1;
    return res;
}

ll exgcd(ll &x, ll &y, ll a, ll b) {
    if (!b) {
        x = 1, y = 0;
        return a;
    }
    ll t = exgcd(y, x, b, a % b);
    y -= x * (a / b);
    return t;
}

ll BSGS(ll a, ll b, ll mod, ll q) {
    H.clear();
    ll Q, p = ceil(sqrt(mod)), x, y;
    exgcd(x, y, q, mod), b = (b * x % mod + mod) % mod,
                         Q = expow(a, p, mod), exgcd(x, y, Q, mod), Q = (x % mod + mod) % mod;
    for (ll i = 1, j = 0; j <= p; ++j, i = i * a % mod) if (!H.count(i)) H[i] = j;
    for (ll i = b, j = 0; j <= p; ++j, i = i * Q % mod) if (H[i]) return j * p + H[i];
    return -1;
}

ll exBSGS(ll N, ll M, ll P) {
    ll q = 1;
    ll k = 0, ret;
    if (M == 1) return 0;
    while ((ret = gcd(N, P)) > 1) {
        if (M % ret) return -1;
        ++k, M /= ret, P /= ret, q = q * (N / ret) % P;
        if (q == M) return k;
    }
    return (ret = BSGS(N, M, P, q)) == -1 ? -1 : ret + k;
}

int main() {
    while (true) {
        int N, M, P;
        scanf("%d%d%d", &N, &P, &M);
        if (!N && !M && !P) break;
        N %= P, M %= P;
        int ans = exBSGS(N, M, P);
        if (ans == -1)
            printf("No Solution\n");
        else
            printf("%d\n", ans);
    }
}