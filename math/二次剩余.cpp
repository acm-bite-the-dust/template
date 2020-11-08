/*
 * 二次剩余，mod为奇素数时有解
 * 解最多两个，为相反数x,(mod-x)
 * n为0特判
 * 算法返回-1则无解
 * mod为2返回n
 */
#include<bits/stdc++.h>

using namespace std;
typedef long long ll;
ll w;
struct num {
    ll x, y;
};

num mul(num a, num b, ll p) {
    num ans = {0, 0};
    ans.x = ((a.x * b.x % p + a.y * b.y % p * w % p) % p + p) % p;
    ans.y = ((a.x * b.y % p + a.y * b.x % p) % p + p) % p;
    return ans;
}

ll qpow_real(ll a, ll b, ll p) {
    ll ans = 1;
    while (b) {
        if (b & 1)ans = 1ll * ans % p * a % p;
        a = a % p * a % p;
        b >>= 1;
    }
    return ans % p;
}

ll qpow_imag(num a, ll b, ll p) {
    num ans = {1, 0};
    while (b) {
        if (b & 1)ans = mul(ans, a, p);
        a = mul(a, a, p);
        b >>= 1;
    }
    return ans.x % p;
}

ll solve(ll n, ll p) {
    n %= p;
    if (p == 2)return n;
    if (qpow_real(n, (p - 1) / 2, p) == p - 1)return -1;//不存在
    ll a;
    while (1) {
        a = rand() % p;
        w = ((a * a % p - n) % p + p) % p;
        if (qpow_real(w, (p - 1) / 2, p) == p - 1)break;
    }
    num x = {a, 1};
    return qpow_imag(x, (p + 1) / 2, p);
}

int main() {
    srand(time(0));
    int t;
    scanf("%d", &t);
    while (t--) {
        ll n, p;
        scanf("%lld%lld", &n, &p);
        if (!n) {
            printf("0\n");
            continue;
        }
        ll ans1 = solve(n, p), ans2;
        if (ans1 == -1) {
            printf("Hola!\n");
        } else {
            ans2 = p - ans1;
            if (ans1 > ans2)swap(ans1, ans2);
            if (ans1 == ans2)printf("%lld\n", ans1);
            else printf("%lld %lld\n", ans1, ans2);
        }
    }
}