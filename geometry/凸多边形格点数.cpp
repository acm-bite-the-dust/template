//计算凸多边形格点数
#include <iostream>
#include <algorithm>

using namespace std;

typedef long long ll;

const int maxn = 1e5;

ll x[maxn + 1] = {};
ll y[maxn + 1] = {};

ll calc(int a, int b) {
    if (y[a] == y[b]) return abs(x[a] - x[b]) - 1;
    if (x[a] == x[b]) return abs(y[a] - y[b]) - 1;
    return __gcd(abs(y[a] - y[b]), abs(x[a] - x[b])) - 1;
}

int main() {
    int n;
    scanf("%d", &n);
    for (int i = 0; i < n; i++) scanf("%lld%lld", &x[i], &y[i]);
    for (int i = 1; i < n; i++) {
        x[i] -= x[0];
        y[i] -= y[0];
    }
    x[n] = y[n] = x[0] = y[0] = 0;
    //2c=2s-a-b+2
    ll ans = 0;
    ll sum = 0;
    for (int i = 1; i <= n; i++) {
        sum += x[i] * y[i + 1] - y[i] * x[i + 1];
    }
    ll b = 0;
    for (int i = 0; i < n; i++) {
        b += calc(i, i + 1);
    }
    cout << (abs(sum) - n - b + 2) / 2;
    return 0;
}