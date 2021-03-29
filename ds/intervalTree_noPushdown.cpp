#include <bits/stdc++.h>
using namespace std;
using ll = long long;
const int MAXN = 100005;
ll tree[MAXN*4];
ll lazy[MAXN*4];
ll arr[MAXN];
void build(int rt, int l, int r) {
    if (l == r) {
        tree[rt] = arr[l];
        return;
    }
    int m = (l + r) / 2;
    build(rt*2, l, m);
    build(rt*2+1, m+1, r);
    tree[rt] = tree[rt*2] + tree[rt*2+1];
}
void update(int rt, int l, int r, int L, int R, ll x) {
    ll len = min(R, r) - max(L, l) + 1;
    tree[rt] += x * len;
    if (L <= l && r <= R) {
        lazy[rt] += x;
        return;
    }
    int m = (l + r) / 2;
    if (L <= m) update(rt*2, l, m, L, R, x);
    if (R > m) update(rt*2+1, m+1, r, L, R, x);
}
ll query(int rt, int l, int r, int L, int R, ll x) {
    ll len = min(R, r) - max(L, l) + 1;
    if (L <= l && r <= R) return tree[rt] + x * len;
    int m = (l + r) / 2;
    ll ans = 0;
    if (L <= m) ans += query(rt*2, l, m, L, R, x+lazy[rt]);
    if (R > m) ans += query(rt*2+1, m+1, r, L, R, x+lazy[rt]);
    return ans;
}
int main() {
//    freopen("../test.in", "r", stdin);
    int n, m; scanf("%d%d", &n, &m);
    for (int i = 1; i <= n; i++) scanf("%lld", arr+i);
    build(1, 1, n);
    while (m--) {
        int opt; scanf("%d", &opt);
        if (opt == 1) {
            int x, y, k; scanf("%d%d%d", &x, &y, &k);
            update(1, 1, n, x, y, k);
        } else {
            int x, y; scanf("%d%d", &x, &y);
            ll ans = query(1, 1, n, x, y, 0);
            printf("%lld\n", ans);
        }
    }
}