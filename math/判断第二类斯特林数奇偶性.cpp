#include <bits/stdc++.h>

using namespace std;
typedef long long ll;

int solve(ll n) {
    int ret = 0;
    while (n) {
        n >>= 1;
        ret += n;
    }
    return ret;
}

int main() {
    ll n, m;
    cin >> n >> m;
    ll a = n - m;
    ll b = (m - 1) / 2;
    ll ans1 = solve(a + b);
    ll ans2 = solve(a);
    ll ans3 = solve(b);
    if (ans1 == ans2 + ans3) {
        printf("1\n");
    } else {
        printf("0\n");
    }
    return 0;
}