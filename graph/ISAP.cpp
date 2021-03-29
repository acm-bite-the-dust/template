//ISAP算法,优化的Dinic,复杂度仍按O(n^2m)计算,大多数情况下比HLPP快(HLPP复杂度：O(n^2sqrt(m)))
#include <bits/stdc++.h>

using namespace std;
const int inf = 0x3f3f3f3f;
const int N = 1210;

//quick read
template<typename T>
inline void read(T &x) {
    int s = 1;
    x = 0;
    char ch = getchar();
    while (ch < '0' || ch > '9') {
        if (ch == '-') s = -1;
        ch = getchar();
    }
    while (ch >= '0' && ch <= '9') {
        x = (x << 3) + (x << 1) + (ch ^ 48);
        ch = getchar();
    }
    x *= s;
}

template<typename T, typename... Args>
inline void read(T &x, Args &... args) {
    read(x);
    read(args...);
}


int tot = 1, head[N];
int n, m, s, t;

struct node {
    int v;
    int next;
    int flow;
} edge[240010];

inline void addedge(int u, int v, int flow) {
    edge[++tot].v = v;
    edge[tot].flow = flow;
    edge[tot].next = head[u];
    head[u] = tot;
}

int dep[N], gap[N];

void bfs() {
    memset(dep, -1, sizeof(dep));
    memset(gap, 0, sizeof(gap));
    dep[t] = 0;
    gap[0] = 1;
    queue<int> q;
    q.push(t);
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        for (int i = head[u]; i; i = edge[i].next) {
            int v = edge[i].v;
            if (dep[v] != -1)continue;
            q.push(v);
            dep[v] = dep[u] + 1;
            gap[dep[v]]++;
        }
    }
}

int maxflow;

int dfs(int u, int flow) {
    if (u == t) {
        maxflow += flow;
        return flow;
    }
    int used = 0;
    for (int i = head[u]; i; i = edge[i].next) {
        int d = edge[i].v;
        if (edge[i].flow && dep[d] + 1 == dep[u]) {
            int mi = dfs(d, min(edge[i].flow, flow - used));
            if (mi) {
                edge[i].flow -= mi;
                edge[i ^ 1].flow += mi;
                used += mi;
            }
            if (used == flow)return used;
        }
    }
    --gap[dep[u]];
    if (gap[dep[u]] == 0)dep[s] = n + 1;
    dep[u]++;
    gap[dep[u]]++;
    return used;
}

int ISAP() {
    maxflow = 0;
    bfs();
    while (dep[s] < n)dfs(s, inf);
    return maxflow;
}

int main() {
    read(n, m, s, t);
    for (int i = 1; i <= m; i++) {
        int u, v, flow;
        read(u, v, flow);
        addedge(u, v, flow);
        addedge(v, u, 0);
    }
    printf("%d", ISAP());
    return 0;
}