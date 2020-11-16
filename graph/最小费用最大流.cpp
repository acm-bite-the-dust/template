//EK算法求解最小费用最大流
#include<iostream>
#include<cstring>
#include<algorithm>
#include<queue>

using namespace std;

const int maxn = 5010;
const int maxe = 50010;
int top = 1, head[maxn];
const int inf = 1 << 30;

int dist[maxn];
int vis[maxn];

struct Node {
    int v;
    int w;//该边单位流量的费用
    int next;
    int val;//该边最大流量
} node[maxe * 2];

struct P {
    int fa;//增广路上该点的父亲
    int edge;//该点与父亲相连的边的编号
} pre[maxn];

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

void init() {
    top = 1;
    memset(head, 0, sizeof(head));
}

inline void addedge(int u, int v, int val, int w) {
    node[++top].v = v;
    node[top].val = val;
    node[top].w = w;
    node[top].next = head[u];
    head[u] = top;
}

inline void add2edge(int u, int v, int val, int w) {
    addedge(u, v, val, w);
    addedge(v, u, 0, -w);
}

bool spfa(int s, int t) {
    memset(pre, 0, sizeof(pre));
    memset(dist, 0x3f, sizeof(dist));
    memset(vis, 0, sizeof(vis));
    queue<int> q;
    q.push(s);
    vis[s] = 1;
    dist[s] = 0;
    while (!q.empty()) {
        int u = q.front();
        vis[u] = 0;
        q.pop();
        int i, d, w;
        for (i = head[u]; i; i = node[i].next) {
            d = node[i].v;
            w = node[i].w;
            if (node[i].val > 0 && dist[d] > dist[u] + w) {
                dist[d] = dist[u] + w;
                pre[d].fa = u;
                pre[d].edge = i;
                if (vis[d] == 0) {
                    q.push(d);
                    vis[d] = 1;
                }
            }
        }
    }
    return dist[t] != 0x3f3f3f3f;
}//寻找费用最小的增广路

auto EK(int s, int t) {//最小费用最大流的EK算法
    int maxflow = 0;
    int cost = 0;
    int mi;
    int i;
    while (spfa(s, t)) {
        mi = inf;
        for (i = t; i != s; i = pre[i].fa)mi = min(mi, node[pre[i].edge].val);
        for (i = t; i != s; i = pre[i].fa) {
            node[pre[i].edge].val -= mi;
            node[pre[i].edge ^ 1].val += mi;
        }
        maxflow += mi;
        cost += mi * dist[t];
    }
    return make_pair(maxflow, cost);
}

int main() {
    init();
    int n, m, s, t;
    read(n, m, s, t);
    for (int i = 0; i < m; i++) {
        int u, v, w, f;
        read(u, v, w, f);
        add2edge(u, v, w, f);
    }
    auto res = EK(s, t);
    printf("%d %d", res.first, res.second);
    return 0;
}