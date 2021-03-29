#include<cstdio>
#include<algorithm>
#include<queue>

using namespace std;

const int inf = 0x3f3f3f3f;

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

const int maxn = 5000 + 5;
const int maxe = 50000 * 2 + 5;

int num = 0;
int head[maxn];
struct node {
    int v, flow, dis;
    int next;
} edge[maxe];

inline void addedge(int x, int y, int w, int dis) {
    edge[num].v = y;
    edge[num].flow = w;
    edge[num].dis = dis;
    edge[num].next = head[x];
    head[x] = num++;
}

inline void add2edge(int x, int y, int w, int dis) {
    addedge(x, y, w, dis);
    addedge(y, x, 0, -dis);
}

int flow[maxn], dis[maxn];
int pre[maxn], last[maxn];
bool vis[maxn];

bool spfa(int n, int s, int t) {
    for (int i = 1; i <= n; ++i) flow[i] = dis[i] = inf, vis[i] = false;
    vis[s] = true;
    dis[s] = 0;
    pre[t] = -1;
    queue<int> q;
    q.push(s);
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        vis[u] = false;
        for (int i = head[u]; i != -1; i = edge[i].next) {
            int v = edge[i].v;
            if (edge[i].flow && dis[v] > dis[u] + edge[i].dis) {
                dis[v] = dis[u] + edge[i].dis;
                pre[v] = u;
                last[v] = i;
                flow[v] = min(flow[u], edge[i].flow);
                if (!vis[v]) {
                    vis[v] = true;
                    q.push(v);
                }
            }
        }
    }
    return pre[t] != -1;
}

pair<int, int> dinic(int n, int s, int t) {
    int maxflow = 0, mincost = 0;
    while (spfa(n, s, t)) {
        int u = t;
        maxflow += flow[u];
        mincost += flow[u] * dis[u];
        while (u != s) {
            edge[last[u]].flow -= flow[t];
            edge[last[u] ^ 1].flow += flow[t];
            u = pre[u];
        }
    }
    return {maxflow, mincost};
}

void init(int n){
	num = 0;
    for (int i = 1; i <= n; ++i) head[i] = -1;
}

int main() {
    int n, m, s, t;
    read(n, m, s, t);
    for (int i = 0; i != m; ++i) {
        int x, y, w, f;
        read(x, y, w, f);
        add2edge(x, y, w, f);
    }
    auto res = dinic(n, s, t);
    printf("%d %d\n", res.first, res.second);
    return 0;
}
