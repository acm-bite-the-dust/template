#include <iostream>
#include <algorithm>

using namespace std;

const int maxn = 500010;
struct node {
    int t, nex;
} e[maxn * 2];

int head[maxn], tot;

void addedge(int x, int y) {
    e[++tot].t = y;
    e[tot].nex = head[x];
    head[x] = tot;
}

int dep[maxn], dp[maxn][22], lg[maxn];

void dfs(int now, int fa) {
    dp[now][0] = fa;
    dep[now] = dep[fa] + 1;
    for (int i = 1; i <= lg[dep[now]]; ++i)
        dp[now][i] = dp[dp[now][i - 1]][i - 1];
    for (int i = head[now]; i; i = e[i].nex)
        if (e[i].t != fa) dfs(e[i].t, now);
}

int LCA(int x, int y) {
    if (dep[x] < dep[y]) swap(x, y);
    while (dep[x] > dep[y])
        x = dp[x][lg[dep[x] - dep[y]] - 1];
    if (x == y) return x;
    for (int k = lg[dep[x]] - 1; k >= 0; --k)
        if (dp[x][k] != dp[y][k])
            x = dp[x][k], y = dp[y][k];
    return dp[x][0];
}

int main() {
    int n, m, s;
    scanf("%d%d%d", &n, &m, &s);
    for (int i = 1; i <= n - 1; ++i) {
        int x, y;
        scanf("%d%d", &x, &y);
        addedge(x, y);
        addedge(y, x);
    }
    for (int i = 1; i <= n; ++i)
        lg[i] = lg[i - 1] + (1 << lg[i - 1] == i);
    dfs(s, 0);
    for (int i = 1; i <= m; ++i) {
        int x, y;
        scanf("%d%d", &x, &y);
        printf("%d\n", LCA(x, y));
    }
    return 0;
}
