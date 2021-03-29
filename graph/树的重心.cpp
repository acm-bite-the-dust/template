#include <cstdio>
#include <vector>
#include <algorithm>

using namespace std;

const int maxn = 1e5 + 10;

int num = 0;
int head[maxn];
struct node {
    int v, next;
} edge[maxn];

inline void addedge(int x, int y) {
    edge[num].v = y;
    edge[num].next = head[x];
    head[x] = num++;
}

int dp[maxn] = {};

int n;

int solve(int x, int fa) {
    int tot = 0;
    for (int i = head[x]; i != -1; i = edge[i].next) {
        int c = edge[i].v;
        if (c != fa) {
            int res = solve(c, x);
            dp[x] = max(dp[x], res);
            tot += res;
        }
    }
    dp[x] = max(dp[x], n - 1 - tot);
    return tot + 1;
}

void init() {
    num = 0;
    for (int i = 1; i <= n; i++) head[i] = -1, dp[i] = 0;
}

int main() {
    while (~scanf("%d", &n)) {
        init();
        for (int i = 1; i < n; i++) {
            int a, b;
            scanf("%d%d", &a, &b);
            addedge(a, b);
            addedge(b, a);
        }
        solve(1, 1);
        int m = dp[1];
        for (int i = 2; i <= n; i++) m = min(m, dp[i]);
        for (int i = 1; i <= n; i++) {
            if (dp[i] == m) printf("%d ", i);
        }
        printf("\n");
    }
    return 0;
}