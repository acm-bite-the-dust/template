//EK算法求解最大流
#include<iostream>
#include<cstring>
#include<queue>

using namespace std;
const int inf = 1 << 30;

const int maxn = 1e5;

struct Node {
    int v;
    int val;
    int next;
} node[2 * maxn + 1];
int top = 1, head[maxn + 1];//top必须从一个奇数开始，一般用-1但我不习惯，解释见下方

void init() {
    top = 1;
    memset(head, 0, sizeof(head));
}

inline void addedge(int u, int v, int val) {
    node[++top].v = v;
    node[top].val = val;
    node[top].next = head[u];
    head[u] = top;
}

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

int vis[maxn + 1];//点是访问过里

struct Pre {
    int v;//该点的前一个点（从起点过来）
    int edge;//与该点相连的边（靠近起点的）
} pre[maxn + 1];

inline bool bfs(int s, int t) {
    queue<int> q;
    memset(vis, 0, sizeof(vis));
    memset(pre, -1, sizeof(pre));
    vis[s] = 1;
    q.push(s);
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        for (int i = head[u]; i; i = node[i].next) {
            int d = node[i].v;
            if (!vis[d] && node[i].val) {//node[i].val==0则已经该路径满了
                pre[d].v = u;
                pre[d].edge = i;
                if (d == t)return true;
                vis[d] = 1;
                q.push(d);
            }
        }
    }
    return false;
}//是否有增广路
int EK(int s, int t) {
    int ans = 0;
    while (bfs(s, t)) {
        int mi = inf;
        for (int i = t; i != s; i = pre[i].v) {
            mi = min(mi, node[pre[i].edge].val);//每次只能增加增广路上最小的边的权值
        }
        for (int i = t; i != s; i = pre[i].v) {
            node[pre[i].edge].val -= mi;
            node[pre[i].edge ^ 1].val += mi;
            //反向的边的编号是正向边的编号^1
            //这就是为什么top开始时必须是奇数
        }
        ans += mi;
    }
    return ans;
}

int main() {
    init();
    int n, m, s, t;
    read(n, m, s, t);
    int u, v, w;
    for (int i = 1; i <= m; i++) {
        read(u, v, w);
        addedge(u, v, w);
        addedge(v, u, 0);
    }
    printf("%d\n", EK(s, t));
    return 0;
}