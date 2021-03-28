// Problem: P3384 【模板】轻重链剖分
// Contest: Luogu
// URL: https://www.luogu.com.cn/problem/P3384
// Memory Limit: 125 MB
// Time Limit: 1000 ms
// 已知一棵包含 N 个结点的树（连通且无环），每个节点上包含一个数值，需要支持以下操作：
// 操作 1： 格式： 1 x y z 表示将树从 x 到 y 结点最短路径上所有节点的值都加上 z
// 操作 2： 格式： 2 x y 表示求树从 x 到 y 结点最短路径上所有节点的值之和
// 操作 3： 格式： 3 x z 表示将以 x 为根节点的子树内所有节点值都加上 z
// 操作 4： 格式： 4 x 表示求以 x 为根节点的子树内所有节点值之和
// Powered by CP Editor (https://cpeditor.org)

#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
const int N = 1e5 + 10;

int n, q, rt;
int dp[N], fa[N], son[N], ord[N], dep[N], a[N], top[N];
ll val[N], mod;
vector<int> G[N];

ll T[N << 2];
ll lazy[N << 2];

void push_down(int l, int r, int id) {
    if (lazy[id] == 0) return;
    int m = (l + r) >> 1;
    T[id << 1] += (m - l + 1) * lazy[id] % mod;
    T[id << 1] %= mod;
    T[id << 1 | 1] += (r - m) * lazy[id] % mod;
    T[id << 1 | 1] %= mod;
    lazy[id << 1] += lazy[id];
    lazy[id << 1] %= mod;
    lazy[id << 1 | 1] += lazy[id];
    lazy[id << 1 | 1] %= mod;
    lazy[id] = 0;
}

void push_up(int id) {
    T[id] = T[id << 1] + T[id << 1 | 1];
    T[id] %= mod;
}

void update(int ql, int qr, int l, int r, ll k, int id) {
    if (ql <= l && r <= qr) {
        lazy[id] += k;
        lazy[id] %= mod;
        T[id] += (r - l + 1) * k % mod;
        T[id] %= mod;
        return;
    }
    push_down(l, r, id);
    int m = (l + r) >> 1;
    if (ql <= m) update(ql, qr, l, m, k, id << 1);
    if (qr > m) update(ql, qr, m + 1, r, k, id << 1 | 1);
    push_up(id);
}

ll query(int ql, int qr, int l, int r, int id) {
    if (ql <= l && r <= qr) {
        return T[id] % mod;
    }
    push_down(l, r, id);
    ll ret = 0;
    int m = (l + r) >> 1;
    if (ql <= m) ret += query(ql, qr, l, m, id << 1) % mod, ret %= mod;
    if (qr > m) ret += query(ql, qr, m + 1, r, id << 1 | 1) % mod, ret %= mod;
    return ret;
}

void build(int l, int r, int id) {
    if (l == r) {
        T[id] = val[l] % mod;
        return;
    }
    int m = (l + r) >> 1;
    build(l, m, id << 1);
    build(m + 1, r, id << 1 | 1);
    push_up(id);
}

void init(int now, int father, int depth) {
    fa[now] = father;
    dp[now] = 1;
    dep[now] = depth;
    int mx = 0;
    for (int i : G[now]) {
        if (i == father) continue;
        init(i, now, depth + 1);
        dp[now] += dp[i];
        mx = max(mx, dp[i]);
    }
    for (int i : G[now]) {
        if (i == father) continue;
        if (dp[i] == mx) {
            son[now] = i;
            break;
        }
    }
}

int cnt;

void dfs(int now, int father, int st) {
    top[now] = st;
    ord[now] = ++cnt;
    val[ord[now]] = a[now];
    if (son[now]) dfs(son[now], now, st);
    for (int i : G[now]) {
        if (i == father) continue;
        if (i == son[now]) continue;
        dfs(i, now, i);
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin >> n >> q >> rt >> mod;
    for (int i = 1; i <= n; ++i) cin >> a[i];
    for (int i = 1; i < n; ++i) {
        int u, v;
        cin >> u >> v;
        G[u].push_back(v);
        G[v].push_back(u);
    }
    init(rt, -1, 1);
    dfs(rt, -1, rt);
    build(1, n, 1);
    while (q--) {
        int op;
        cin >> op;
        if (op == 1) {
            int x, y, z;
            cin >> x >> y >> z;
            z %= mod;
            while (top[y] != top[x]) {
                if (dep[top[x]] > dep[top[y]]) swap(x, y);
                update(ord[top[y]], ord[y], 1, n, z, 1);
                y = fa[top[y]];
            }
            if (dep[x] > dep[y]) swap(x, y);
            update(ord[x], ord[y], 1, n, z, 1);
        } else if (op == 2) {
            int x, y;
            cin >> x >> y;
            ll ans = 0;
            while (top[y] != top[x]) {
                if (dep[top[x]] > dep[top[y]]) swap(x, y);
                ans += query(ord[top[y]], ord[y], 1, n, 1) % mod;
                ans %= mod;
                y = fa[top[y]];
            }
            if (dep[x] > dep[y]) swap(x, y);
            ans += query(ord[x], ord[y], 1, n, 1);
            ans %= mod;
            cout << ans << '\n';
        } else if (op == 3) {
            int x, z;
            cin >> x >> z;
            z %= mod;
            update(ord[x], ord[x] + dp[x] - 1, 1, n, z, 1);
        } else {
            int x;
            cin >> x;
            ll ans = query(ord[x], ord[x] + dp[x] - 1, 1, n, 1) % mod;
            cout << ans << '\n';
        }
    }
    return 0;
}