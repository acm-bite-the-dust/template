//点分治
#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
const int maxn = 20005;
const int inf = 0x3f3f3f3f;
const int mod = 1e9 + 7;

struct edge {
    int to, val;
};
vector<edge> mp[maxn];

int mini, rt, totSZ;
int sz[maxn], dis[maxn];
bool vis[maxn];

int l, r, q[maxn];  //q为每次得到的距离合集

int n, ans; //ans记录合法点对

void getRT(int u, int pre) {    //每次调用 getRT() 前使 mini=inf, totSZ = sz[v];
    sz[u] = 1;
    int mxSub = 0;
    for (auto it: mp[u]) {
        int v = it.to;
        if (v == pre || vis[v]) continue;

        getRT(v, u);
        sz[u] += sz[v];
        mxSub = max(mxSub, sz[v]);
    }

    int mx = max(mxSub, totSZ - sz[u]);
    if (mx < mini) {
        mini = mx;
        rt = u;
    }
}

void getDIS(int u, int pre) {
    q[++r] = dis[u];
    for (auto it: mp[u]) {
        int v = it.to, val = it.val;
        if (v == pre || vis[v]) continue;
        dis[v] = dis[u] + val;
        getDIS(v, u);
    }
}

int calc(int u, int val) {
    l = 1, r = 0;
    dis[u] = val;
    getDIS(u, 0);

    //按照题意处理q

    return sum;
}

void dfs(int u) {
    vis[u] = true;
    ans += calc(u, 0);
    for (auto it: mp[u]) {
        int v = it.to, val = it.val;
        if (vis[v]) continue;

        ans -= calc(v, val);

        mini = inf;
        totSZ = sz[v];
        getRT(v, 0);
        dfs(rt);
    }
}

int main() {

    while (~scanf("%d", &n)) {
        for (int i = 1; i <= n; i++) {
            mp[i].clear();
            vis[i] = false;
        }
        ans = 0;

        for (int i = 1; i < n; i++) {
            int u, v, val;
            scanf("%d%d%d", &u, &v, &val);
            mp[u].push_back(edge{v, val});
            mp[v].push_back(edge{u, val});
        }

        mini = inf;
        totSZ = n;
        getRT(1, 0);
        dfs(rt);

        printf("%d\n", ans);
    }

    return 0;
}