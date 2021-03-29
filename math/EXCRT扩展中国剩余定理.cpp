//x % A[i] = B[i] O(nlogn) x为最小解 x + k * lcm都可行
#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
const int N = 1e5 + 10;

ll mul(ll a, ll b, ll mod) {
    ll ret = 0;
    while (b) {
        if (b & 1) ret = (ret + a) % mod;
        a = (a + a) % mod;
        b >>= 1;
    }
    return ret;
}

ll exgcd(ll a, ll b, ll &x, ll &y) {
    ll ret, tmp;
    if (!b) {
        x = 1;
        y = 0;
        return a;
    }
    ret = exgcd(b, a % b, x, y);
    tmp = x;
    x = y;
    y = tmp - a / b * y;
    return ret;
}

ll A[N], B[N];

ll excrt(int n) {
    ll x, y;
    ll M = A[1], ans = B[1];
    for (int i = 2; i <= n; ++i) {
        ll a = M, b = A[i], c = (B[i] - ans % b + b) % b;
        ll g = exgcd(a, b, x, y), bg = b / g;
        if (c % g) return -1;
        x = mul(x, c / g, bg); //可能溢出
        ans += x * M;
        M *= bg;
        ans = (ans % M + M) % M;
    }
    return (ans % M + M) % M;
}

int main() {
    int n;
    cin >> n;
    for (int i = 1; i <= n; ++i) {
        cin >> A[i] >> B[i];
    }
    cout << excrt(n);
    return 0;
}