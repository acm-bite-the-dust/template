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

const int maxn = 1e5 + 5;
const int maxe = 1e6;

int num = 0;
int head[maxn];
struct node {
    int v, flow;
    int next;
} edge[2 * maxe + 10];

inline void addedge(int x, int y, int w) {
    edge[num].v = y;
    edge[num].flow = w;
    edge[num].next = head[x];
    head[x] = num++;
}

inline void add2edge(int x, int y, int w) {
    addedge(x, y, w);
    addedge(y, x, 0);
}

int cur[maxn], dis[maxn];

bool bfs(int n, int s, int t) {
    for (int i = 1; i <= n; ++i) dis[i] = 0;
    dis[s] = 1;
    queue<int> q;
    q.push(s);
    while (!q.empty()) {
        int u = q.front();
        if (u == t) return true;
        q.pop();
        for (int i = head[u]; i != -1; i = edge[i].next) {
            int v = edge[i].v;
            if (!dis[v] && edge[i].flow) {
                dis[v] = dis[u] + 1;
                q.push(v);
            }
        }
    }
    return dis[t] != 0;
}

int dfs(int u, int f, int t) { //��ǰ�ڵ㣬��ĿǰΪֹ���������������ؿ���������
    if (u == t || (!f)) return f; // ������ || ���޲���
    int flow = 0;
    for (int &i = cur[u]; i != -1; i = edge[i].next) {// '&'Ŀ�� -> ��ǰ���Ż�
        int v = edge[i].v;
        if (dis[v] == dis[u] + 1) {
            int di = dfs(v, min(f, edge[i].flow), t);
            if (di > 0) {
                flow += di;//���ﲻ����return�����Ǽ�¼һ������������������������ٽ�������һ����
                f -= di;//�Ѵ�u��ʼ�������������Ӧ��ȥ
                edge[i].flow -= di;
                edge[i ^ 1].flow += di;
                if (!f) break;// û������
            }
        }
    }
    if (!flow) dis[u] = -1;//���������û�������һ������,ը����
    return flow;
}

int dinic(int n, int s, int t) {
    int ans = 0;
    while (bfs(n, s, t)) {
        for (int i = 1; i <= n; ++i) cur[i] = head[i];
        ans += dfs(s, inf, t);
    }
    return ans;
}

void init(int n) {
	num = 0;
	for (int i = 1; i <= n; ++i) head[i] = -1;
}

int main() {
    int n, m, s, t;
    read(n, m, s, t);
    init(n);
    for (int i = 0; i != m; ++i) {
        int x, y, w;
        read(x, y, w);
        add2edge(x, y, w);
    }
    printf("%d\n", dinic(n, s, t));
    return 0;
}
