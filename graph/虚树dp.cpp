// Problem: P2495 [SDOI2011]消耗战
// 1为根，每次可以炸毁一个边，花费为边权。
// m次询问，每次要求让k个关键点无法到达1，问最小花费
// 对于每次询问单独建树，dfs序+单调栈优化

#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
const int N = 3e5 + 10;
vector<pair<int, ll>> G[N];
vector<int> g[N];
int ord[N]; //dfs序
int stk[N];
int dep[N];
int parent[N][30];
bool inq[N]; //是否为关键点
ll minV[N]; //从1到i不联通的最小花费
ll dp[N]; //ans
int now;
int maxx = 0;

void dfs1(int u, int fa) {
    ++now;
    ord[u] = now;
    dep[u] = dep[fa] + 1;
    maxx = max(maxx, dep[u]);
    parent[u][0] = fa;
    for (auto pii : G[u]) {
        int v = pii.first;
        ll w = pii.second;
        if (v == fa) continue;
        minV[v] = min(minV[u], w);
        dfs1(v, u);
    }
}

void dfs2(int u) {
    ll sum = 0;
    for (int v : g[u]) {
        dfs2(v);
        sum += dp[v];
    }
    if (inq[u]) {
        dp[u] = minV[u];
    } else {
        dp[u] = min(minV[u], sum);
    }
    inq[u] = false;
    g[u].clear();
}

int n;

void init() {
    memset(minV, 0x3f, sizeof minV);
    dfs1(1, 0);
    for (int i = 1; (1 << i) <= n; ++i) {
        for (int j = 1; j <= n; ++j) {
            if (parent[j][i - 1] == 0) parent[j][i] = 0;
            else parent[j][i] = parent[parent[j][i - 1]][i - 1];
        }
    }
}

int LCA(int ta, int tb) {
    if (dep[ta] < dep[tb]) swap(ta, tb);
    for (int i = 0; dep[ta] != dep[tb]; ++i) {
        if ((dep[ta] - dep[tb]) >> i & 1)
            ta = parent[ta][i];
    }
    if (ta == tb) return ta;
    for (int i = log2(maxx); i >= 0; --i) {
        if (parent[ta][i] != parent[tb][i]) {
            ta = parent[ta][i];
            tb = parent[tb][i];
        }
    }
    return parent[ta][0];
}

void build() {
    vector<pair<int, int>> v;
    int k;
    cin >> k;
    for (int i = 0; i < k; ++i) {
        int id;
        cin >> id;
        inq[id] = true;
        v.emplace_back(ord[id], id);
    }
    sort(v.begin(), v.end());
    now = 0;
    stk[++now] = 1;
    g[1].clear();
    for (int i = 0; i < k; ++i) {
        int id = v[i].second;
        if (id != 1) {
            int lca = LCA(id, stk[now]);
            if (lca != stk[now]) {
                while (ord[lca] < ord[stk[now - 1]]) {
                    g[stk[now - 1]].push_back(stk[now]);
                    now--;
                }
                if (ord[lca] > ord[stk[now - 1]]) {
                    g[lca].clear();
                    g[lca].push_back(stk[now]);
                    now--;
                    stk[++now] = lca;
                } else {
                    g[lca].push_back(stk[now]);
                    now--;
                }
            }
        }
        g[id].clear();
        stk[++now] = id;
    }
    for (int i = 1; i < now; ++i) {
        g[stk[i]].push_back(stk[i + 1]);
    }
}

void solve() {
    cin >> n;
    for (int i = 1; i < n; ++i) {
        int u, v, w;
        cin >> u >> v >> w;
        G[u].emplace_back(v, w);
        G[v].emplace_back(u, w);
    }
    init();
    int m;
    cin >> m;
    while (m--) {
        build();
        dfs2(stk[1]);
        printf("%lld\n", dp[stk[1]]);
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    int T = 1;
    // cin >> T;
    while (T--) {
        solve();
    }
    return 0;
}