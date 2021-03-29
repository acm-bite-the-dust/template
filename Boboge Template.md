# Boboge Template

0. 比赛经验谈
	- 做题策略
	- 心态
1. dp
	- 背包
	- 数位
2. ds
	- 线段树
	- 主席树
3. geo
	- 凸多边形格点数
	- 凸包
	- Geo常用&任意多边形面积
4. graph
	- dijkstra
	- ISAP
	- 倍增LCA
	- Tarjan缩点
	- 点分治
	- 树重心
	- 拓扑
	- dinic
5. math
	- BSGS原根
	- EXCRT(含EXGCD)
	- 倍增求矩阵幂和
	- 高斯消元
	- 格点数
	- Miller-Rabin & Pollard Rho
	- 矩阵快速幂
	- 快速判断能否被某个数整除
	- 欧拉降幂
	- 判断第二类斯特林数奇偶性
	- 三分
	- 线性基
6. string
	- 字典树
	- KMP
	- ac自动机
	- 马拉车
	- 回文自动机
7. others
	- 快读
	- Linux对拍脚本
	- Windows对拍脚本



## 0.比赛经验谈

### 比赛策略

先做简单题，后做难题。先做思维，后做模拟。

多花时间在理解、思考、优化，对手速有一定自信。

跟榜跟榜，过的人多的题肯定是有切入点，找切入点。

打自己熟悉的代码，如果是以前完全没见过的类型且没有清晰思路的话，重新审题。

### 心态

卡题的时候先看榜，尽量把题都读完。

怎么说呢，心态这一块必须拿捏住，咱不缺比赛打。

一直以来积累的东西并非全部木大。



### 1.dp

### 背包

![image-20201025153021578](C:\Users\Boboge\AppData\Roaming\Typora\typora-user-images\image-20201025153021578.png)

![image-20201025153039657](C:\Users\Boboge\AppData\Roaming\Typora\typora-user-images\image-20201025153039657.png)

![image-20201025153057665](C:\Users\Boboge\AppData\Roaming\Typora\typora-user-images\image-20201025153057665.png)



![image-20201025153117106](C:\Users\Boboge\AppData\Roaming\Typora\typora-user-images\image-20201025153117106.png)



![image-20201025153128327](C:\Users\Boboge\AppData\Roaming\Typora\typora-user-images\image-20201025153128327.png)



![image-20201025153140547](C:\Users\Boboge\AppData\Roaming\Typora\typora-user-images\image-20201025153140547.png)



![image-20201025153156341](C:\Users\Boboge\AppData\Roaming\Typora\typora-user-images\image-20201025153156341.png)



![image-20201025153215685](C:\Users\Boboge\AppData\Roaming\Typora\typora-user-images\image-20201025153215685.png)



![image-20201025153234184](C:\Users\Boboge\AppData\Roaming\Typora\typora-user-images\image-20201025153234184.png)

![image-20201025153309767](C:\Users\Boboge\AppData\Roaming\Typora\typora-user-images\image-20201025153309767.png)



### 数位dp

```cpp
// Problem: windy 数 相邻数位差至少为2的个数，不含前导零，[l,r]个数

#include <bits/stdc++.h>

#define pb push_back
#define eb emplace_back
#define size(a) (int)((a).size())

using namespace std;
typedef long long ll;
typedef long double ld;
typedef unsigned long long ull;
const int N = 2e5 + 10;
int dp[10][10][2][2];  // len dig lim st
int num[10];

int calc(int x) {
  if (!x) return 0;
  memset(dp, 0, sizeof(dp));
  int tmp = x;
  int len = 0;
  for (int i = 0; i < 10; ++i) {
    num[i] = tmp % 10;
    tmp /= 10;
    if (num[i] != 0) len = i;
  }
  for (int i = 1; i < num[len]; ++i) dp[len][i][0][1] = 1;
  dp[len][num[len]][1][1] = 1;
  for (int i = len - 1; i >= 0; --i) {
    for (int j = 1; j <= 9; ++j) {
      dp[i][j][0][1] = 1;
    }
    for (int j = 0; j <= 9; ++j) {
      for (int k = 0; k <= 9; ++k) {
        if (abs(j - k) < 2) continue;
        if (k == num[i]) {
          dp[i][k][1][0] += dp[i + 1][j][1][1];
          dp[i][k][1][0] += dp[i + 1][j][1][0];
        } else if (k < num[i]) {
          dp[i][k][0][0] += dp[i + 1][j][1][1];
          dp[i][k][0][0] += dp[i + 1][j][1][0];
        }
        dp[i][k][0][0] += dp[i + 1][j][0][1];
        dp[i][k][0][0] += dp[i + 1][j][0][0];
      }
    }
  }
  int ans = 0;
  for (int i = 0; i <= 9; ++i) {
    for (int j = 0; j <= 1; ++j) {
      for (int k = 0; k <= 1; ++k) {
        ans += dp[0][i][j][k];
      }
    }
  }
  return ans;
}

void solve() {
  int a, b;
  cin >> a >> b;
  int ans = calc(b) - calc(a - 1);
  cout << ans << '\n';
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
```



## 2. ds

### 线段树

```cpp
#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
const int N = 1e5 + 10;
ll T[N << 2];
ll lazy[N << 2];

void push_down(int l, int r, int id) {
    int m = (l + r) >> 1;
    T[id << 1] += (m - l + 1) * lazy[id];
    T[id << 1 | 1] += (r - m) * lazy[id];
    lazy[id << 1] += lazy[id];
    lazy[id << 1 | 1] += lazy[id];
    lazy[id] = 0;
}

void update(int ql, int qr, int l, int r, ll k, int id) {
    if (ql <= l && r <= qr) {
        lazy[id] += k;
        T[id] += (r - l + 1) * k;
        return;
    }
    push_down(l, r, id);
    int m = (l + r) >> 1;
    if (ql <= m) update(ql, qr, l, m, k, id << 1);
    if (qr > m) update(ql, qr, m + 1, r, k, id << 1 | 1);
    T[id] = T[id << 1] + T[id << 1 | 1];
}

ll query(int ql, int qr, int l, int r, int id) {
    if (ql <= l && r <= qr) {
        return T[id];
    }
    push_down(l, r, id);
    ll ret = 0;
    int m = (l + r) >> 1;
    if (ql <= m) ret += query(ql, qr, l, m, id << 1);
    if (qr > m) ret += query(ql, qr, m + 1, r, id << 1 | 1);
    return ret;
}

void build(int l, int r, int id) {
    if (l == r) {
        scanf("%d", &T[id]);
        return;
    }
    int m = (l + r) >> 1;
    build(l, m, id << 1);
    build(m + 1, r, id << 1 | 1);
    T[id] = T[id << 1] + T[id << 1 | 1];
}

int main() {
    int n, m;
    scanf("%d%d", &n, &m);
    build(1, n, 1);
    while (m--) {
        int op;
        scanf("%d", &op);
        if (op == 1) {
            int x, y, k;
            scanf("%d%d%d", &x, &y, &k);
            update(x, y, 1, n, k, 1);
        } else {
            int x, y;
            scanf("%d%d", &x, &y);
            printf("%lld\n", query(x, y, 1, n, 1));
        }
    }
    return 0;
}

```



### 主席树

```cpp
//求区间k大
#include <bits/stdc++.h>

using namespace std;
const int maxn = 2e5 + 10;
int tot;
int sum[maxn << 5], id[maxn], ls[maxn << 5], rs[maxn << 5];
int a[maxn], len;


int build(int l, int r) {
    int root = ++tot;
    if (l == r) return root;
    int mid = l + r >> 1;
    ls[root] = build(l, mid);
    rs[root] = build(mid + 1, r);
    return root;
}

int update(int k, int l, int r, int root) {
    int dir = ++tot;
    ls[dir] = ls[root], rs[dir] = rs[root], sum[dir] = sum[root] + 1;
    if (l == r) return dir;
    int mid = l + r >> 1;
    if (k <= mid)
        ls[dir] = update(k, l, mid, ls[dir]);
    else
        rs[dir] = update(k, mid + 1, r, rs[dir]);
    return dir;
}

int query(int u, int v, int l, int r, int k) {
    int mid = l + r >> 1, x = sum[ls[v]] - sum[ls[u]];
    if (l == r) return l;
    if (k <= x)
        return query(ls[u], ls[v], l, mid, k);
    else
        return query(rs[u], rs[v], mid + 1, r, k - x);
}

vector<int> ind;

int getid(int val) {
    return lower_bound(ind.begin(), ind.end(), val) - ind.begin() + 1;
}

int main() {
    int n, m;
    scanf("%d%d", &n, &m);
    for (int i = 1; i <= n; ++i) scanf("%d", a + i), ind.push_back(a[i]);
    sort(ind.begin(), ind.end());
    ind.erase(unique(ind.begin(), ind.end()), ind.end());
    len = ind.size();
    id[0] = build(1, len);
    for (int i = 1; i <= n; ++i) id[i] = update(getid(a[i]), 1, len, id[i - 1]);
    while (m--) {
        int l, r, k;
        scanf("%d%d%d", &l, &r, &k);
        printf("%d\n", ind[query(id[l - 1], id[r], 1, len, k) - 1]);  // 回答询问
    }
    return 0;
}
```



## 3. geo

### 凸多边形格点数

```cpp
//计算凸多边形格点数
#include <iostream>
#include <algorithm>

using namespace std;

typedef long long ll;

const int maxn = 1e5;

ll x[maxn + 1] = {};
ll y[maxn + 1] = {};

ll calc(int a, int b) {
    if (y[a] == y[b]) return abs(x[a] - x[b]) - 1;
    if (x[a] == x[b]) return abs(y[a] - y[b]) - 1;
    return __gcd(abs(y[a] - y[b]), abs(x[a] - x[b])) - 1;
}

int main() {
    int n;
    scanf("%d", &n);
    for (int i = 0; i < n; i++) scanf("%lld%lld", &x[i], &y[i]);
    for (int i = 1; i < n; i++) {
        x[i] -= x[0];
        y[i] -= y[0];
    }
    x[n] = y[n] = x[0] = y[0] = 0;
    //2c=2s-a-b+2
    ll ans = 0;
    ll sum = 0;
    for (int i = 1; i <= n; i++) {
        sum += x[i] * y[i + 1] - y[i] * x[i + 1];
    }
    ll b = 0;
    for (int i = 0; i < n; i++) {
        b += calc(i, i + 1);
    }
    cout << (abs(sum) - n - b + 2) / 2;
    return 0;
}
```



### 凸包

```cpp
//凸包
#include<bits/stdc++.h>
using namespace std;

struct point{
    double x,y;
    point friend operator -(point a,point b)
    {return {a.x-b.x,a.y-b.y};}
}p[105],s[105];
double dis(point a,point b)
{
    point c=a-b;
    return sqrt(c.x*c.x+c.y*c.y);
}
double X(point a,point b)
{
    return a.x*b.y-a.y*b.x;
}
int cmp(point a,point b)
{
    double x=X(a-p[1],b-p[1]);
    
    if(x>0) return 1;
    if(x==0&&dis(a,p[1])<dis(b,p[1])) return 1;
    return 0;
}
double multi(point p1,point p2,point p3)
{
    return X(p2-p1,p3-p1);
}
int main()
{
    int N;
    while(scanf("%d",&N),N)
    {
        for(int i=1;i<=N;i++) cin>>p[i].x>>p[i].y;
        
        if(N==1)
        {
            printf("0.00\n");
            continue;
        }
        else if(N==2)
        {
            printf("%.2lf\n",dis(p[1],p[2]));
            continue;
        }
        
        int k=1;
        for(int i=2;i<=N;i++)
        if(p[i].y<p[k].y||(p[i].y==p[k].y&&p[i].x<p[k].x))k=i;
        swap(p[1],p[k]);
        
        sort(p+2,p+1+N,cmp);
        
        s[1]=p[1];
        s[2]=p[2];
        int t=2;
        for(int i=3;i<=N;i++)
        {
            while(t>=2&&multi(s[t-1],s[t],p[i])<=0) t--;
            s[++t]=p[i];
        }
        double sum=0;
        for(int i=1;i<t;i++)
        {
            sum+=dis(s[i],s[i+1]);
        }
        printf("%.2lf\n",sum+dis(s[1],s[t]));
    }
    return 0;
}
```



### Geo常用&任意多边形面积

```cpp
//geo_base
#include<bits/stdc++.h>
using namespace std;

struct Point
{
	double x,y;
	Point(double x,double y):x(x),y(y){}
	Point operator +(const Point B)const{ return Point(x+B.x,y+B.y); }
	Point operator -(const Point B)const{ return Point(x-B.x,y-B.y); }
	Point operator *(const double B)const{ return Point(x*B,y*B); }
	double Cross(const Point B)const{ return x*B.y-y*B.x; }
	double Dot(const Point B)const{ return x*B.x+y*B.y; }
	double Len()const{ return sqrt(Dot(*this)); }
	double Angle(const Point B)const{ return acos(Dot(B)/Len()/B.Len()); }
};

int dcmp(double x, double eps)
{ //三态函数，克服浮点数精度陷阱，判断x==0?x<0?x>0?
    if (fabs(x) < eps)
        return 0;
    else
        return x < 0 ? -1 : 1;
}

bool SegmentProperIntersection(Point a1, Point a2, Point b1, Point b2)
{
    double c1 = Cross(a2 - a1, b1 - a1), c2 = Cross(a2 - a1, b2 - a1),
           c3 = Cross(b2 - b1, a1 - b1), c4 = Cross(b2 - b1, a2 - b1);
    return dcmp(c1) * dcmp(c2) < 0 && dcmp(c3) * dcmp(c4) < 0;
}

Point GetLineIntersection(Point P, Vector v, Point Q, Vector w)
{
    Vector u = P - Q;
    double t = Cross(w, u) / Cross(v, w);
    return P + v * t;
}

//计算任意多边形的面积，顶点按照顺时针或者逆时针方向排列
double polygon_area(Point *p, int n)
{
    if(n < 3) return 0;

    double sum = 0;
    p[n + 1] = p[1];
    for(int i = 1; i <= n; i++)
        sum += p[i].x * p[i + 1].y - p[i].y * p[i + 1].x;//可以理解为不管这个多边形在哪，都以原点为分割点，就算原点在外面也可以算出，因为有正负可以抵消掉多余的
    sum = fabs(sum / 2.0);
    return sum;
}
```



## 4. Graph

### dijstra

```cpp
//dij
#include <queue>
#include <cstdio>
#include <cstring>
#include <algorithm>
using namespace std;
const int maxn=1e2+7;
const int maxm=1e4+7;
const int inf=1e9+7;
struct edge{
    int to,next,w;
}e[maxm*2];
struct node{
    int id,v;
    bool operator<(const node &t)const{
        if(v!=t.v) return v>t.v;
    }
};
int n,m,dis[maxn],head[maxn],vis[maxn],cnt;
void add(int u,int v,int w){
    e[cnt].next=head[u];
    e[cnt].to=v;
    e[cnt].w=w;
    head[u]=cnt++;
}
void dijkstra(int s){
    for(int i=1;i<=n;i++) dis[i]=inf,vis[i]=0;
    dis[s]=0;
    priority_queue<node>q;
    q.push(node{s,0});
    while(!q.empty()){
        node p=q.top();
        q.pop();
        if(vis[p.id]) continue;
        vis[p.id]=1;
        for(int i=head[p.id];~i;i=e[i].next){
            int v=e[i].to,w=e[i].w;
            if(dis[v]>dis[p.id]+w){
                dis[v]=dis[p.id]+w;
                q.push(node{v,dis[v]});
            }
        }
    }
}
void solve(){
    cnt=0;
    for(int i=1;i<=n;i++) head[i]=-1;
    for(int i=0;i<m;i++){
        int u,v,w;scanf("%d%d%d",&u,&v,&w);
        add(u,v,w);
        add(v,u,w);
    }
    dijkstra(1);
    printf("%d\n",dis[n]);
}
int main(void){
    while(~scanf("%d%d",&n,&m)){
        if(!n&&!m) break;
        solve();
    }
    return 0;
}
```



### ISAP

```cpp
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
```



### 倍增LCA

![LCA 树上倍增](C:\Users\Boboge\Documents\GitHub\template\graph\LCA 树上倍增.png)



### Tarjan缩点

```cpp
const int N = 1e5 + 10;
vector<int> G[N];
int dfn[N], low[N];//dfn[u] -> u被搜索的次序 low[u] -> u子树中dfn最小值(子树包括自身)
int index;//搜索序号
stack<int> S;
bool ins[N];//是否进栈
int col[N], num_color;//染色

void Tarjan(int u) {
    dfn[u] = low[u] = ++index;
    S.push(u);//进栈
    ins[u] = true;
    for (int i :G[u]) {
        int v = i;
        if (!dfn[v]) { //未被访问过
            Tarjan(v);
            low[u] = min(low[u], low[v]); //找爸爸（环开头）最小的
        } else if (ins[v]) { //已被访问过且在栈内,则需要处理; 若不在栈内说明对应强连通分量处理完成
            low[u] = min(low[u], dfn[v]); //判断谁是爸爸
        }
    }
    if (dfn[u] == low[u]) { //发现更新完一轮自己是爸爸(某强连通分量中仅第一个被访问的结点满足dfn[u] == low[u])
        num_color++;
        int tmp;
        do {
            tmp = S.top();
            col[tmp] = num_color; //出栈，染色
            ins[tmp] = false;
            S.pop();
        } while (tmp != u);
    }
}

/*
 * 染好色后, 同一个color可以缩成一个点, 整张图相当于一张DAG
 * DAG可以拓扑, 搜索啥的, 为所欲为
 */
```



### 点分治

```cpp
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
```



### 树的重心

```cpp
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
```



### 拓扑

```cpp
void top() {
    for (int i = 1; i <= n; i++) {
        for (int v : G[i]) {
            in[v]++;
        }
    }
    queue<int> q;
    for (int i = 1; i <= n; i++) {
        if (!in[i]) {
            q.push(i);
        }
    }
    int tot = 0;
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        ans[++tot] = u;
        for (int v : G[u]) {
            in[v]--;
            if (!in[v]) {
                q.push(v);
            }
        }
    }
    if (tot != n) {
        printf("NO\n");
    } else {
        printf("YES\n");
        for (int i = 1; i <= n; i++) {
            printf("%d%c", ans[i], i == n ? '\n' : ' ');
        }
    }
}
```



### Dinic

```cpp
#include <bits/stdc++.h>

using namespace std;
const int N = 5e3 + 10;
const int inf = 0x3f3f3f3f;
int top = 1, head[N];
int n, m, s, t;
int cost, maxflow;
int vis[N];//是否到达过该点 
int dist[N];//到t的单位费用 
struct node {
    int val;
    int v;
    int next;
    int w;
} edge[100010];

void addedge(int u, int v, int val, int w) {
    edge[++top].v = v;
    edge[top].w = w;
    edge[top].val = val;
    edge[top].next = head[u];
    head[u] = top;
}

bool spfa() {
    memset(vis, 0, sizeof(vis));
    memset(dist, 0x3f, sizeof(dist));
    dist[s] = 0;
    vis[s] = 1;
    queue<int> q;
    q.push(s);
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        vis[u] = 0;
        for (int i = head[u]; i; i = edge[i].next) {
            int d = edge[i].v;
            if (dist[d] > dist[u] + edge[i].w && edge[i].val) {
                dist[d] = dist[u] + edge[i].w;
                if (vis[d] == 0) {
                    q.push(d);
                    vis[d] = 1;
                }
            }
        }
    }
    return dist[t] != 0x3f3f3f3f;
}

int dfs(int u, int flow) {
    if (u == t) {
        vis[t] = 1;
        maxflow += flow;
        return flow;
    }//可以到达t，加流 
    int used = 0;//该条路径可用流量 
    vis[u] = 1;
    for (int i = head[u]; i; i = edge[i].next) {
        int d = edge[i].v;
        if ((vis[d] == 0 || d == t) && edge[i].val != 0 && dist[d] == dist[u] + edge[i].w) {
            int minflow = dfs(d, min(flow - used, edge[i].val));
            if (minflow != 0)
                cost += edge[i].w * minflow, edge[i].val -= minflow, edge[i ^ 1].val += minflow, used += minflow;
            //可以到达t，加费用，扣流量 
            if (used == flow)break;
        }
    }
    return used;
}

int mincostmaxflow() {
    while (spfa()) {
        vis[t] = 1;
        while (vis[t]) {
            memset(vis, 0, sizeof(vis));
            dfs(s, inf);
        }
    }
    return maxflow;
}

int main() {
    scanf("%d%d%d%d", &n, &m, &s, &t);
    int u, v, w, val;
    for (int i = 1; i <= m; i++) {
        scanf("%d%d%d%d", &u, &v, &w, &val);
        addedge(u, v, val, w);
        addedge(v, u, 0, -w);
    }
    mincostmaxflow();
    printf("%d %d", maxflow, cost);
    return 0;
}
```



## 5. Math

### BSGS

```cpp
/* exBSGS算法
 * N ^ x = M (mod P), 其中N, P互质
 * 返回的x为最小解
 */
#include<bits/stdc++.h>

using namespace std;
typedef long long ll;
unordered_map<ll, int> H;

ll gcd(ll a, ll b) {
    if (!b) return a;
    return gcd(b, a % b);
}

ll expow(ll a, ll b, ll mod) {
    ll res = 1;
    while (b) res = ((b & 1) ? res * a % mod : res), a = a * a % mod, b >>= 1;
    return res;
}

ll exgcd(ll &x, ll &y, ll a, ll b) {
    if (!b) {
        x = 1, y = 0;
        return a;
    }
    ll t = exgcd(y, x, b, a % b);
    y -= x * (a / b);
    return t;
}

ll BSGS(ll a, ll b, ll mod, ll q) {
    H.clear();
    ll Q, p = ceil(sqrt(mod)), x, y;
    exgcd(x, y, q, mod), b = (b * x % mod + mod) % mod,
                         Q = expow(a, p, mod), exgcd(x, y, Q, mod), Q = (x % mod + mod) % mod;
    for (ll i = 1, j = 0; j <= p; ++j, i = i * a % mod) if (!H.count(i)) H[i] = j;
    for (ll i = b, j = 0; j <= p; ++j, i = i * Q % mod) if (H[i]) return j * p + H[i];
    return -1;
}

ll exBSGS(ll N, ll M, ll P) {
    ll q = 1;
    ll k = 0, ret;
    if (M == 1) return 0;
    while ((ret = gcd(N, P)) > 1) {
        if (M % ret) return -1;
        ++k, M /= ret, P /= ret, q = q * (N / ret) % P;
        if (q == M) return k;
    }
    return (ret = BSGS(N, M, P, q)) == -1 ? -1 : ret + k;
}

int main() {
    while (true) {
        int N, M, P;
        scanf("%d%d%d", &N, &P, &M);
        if (!N && !M && !P) break;
        N %= P, M %= P;
        int ans = exBSGS(N, M, P);
        if (ans == -1)
            printf("No Solution\n");
        else
            printf("%d\n", ans);
    }
}
```



### EXCRT(含EXGCD)

```cpp
//x % A[i] = B[i] O(nlogn) x为最小解 x + k * lcm都可行
#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
const int N = 1e5 + 10;

ll mul(ll a, ll b, ll mod) {
    ll ret = 0;
    while (b) {
        if (b & 1) ret = (ret + a) % mod;
        a = (a + a) % mod;
        b >>= 1;
    }
    return ret;
}

ll exgcd(ll a, ll b, ll &x, ll &y) {
    ll ret, tmp;
    if (!b) {
        x = 1;
        y = 0;
        return a;
    }
    ret = exgcd(b, a % b, x, y);
    tmp = x;
    x = y;
    y = tmp - a / b * y;
    return ret;
}

ll A[N], B[N];

ll excrt(int n) {
    ll x, y;
    ll M = A[1], ans = B[1];
    for (int i = 2; i <= n; ++i) {
        ll a = M, b = A[i], c = (B[i] - ans % b + b) % b;
        ll g = exgcd(a, b, x, y), bg = b / g;
        if (c % g) return -1;
        x = mul(x, c / g, bg); //可能溢出
        ans += x * M;
        M *= bg;
        ans = (ans % M + M) % M;
    }
    return (ans % M + M) % M;
}

int main() {
    int n;
    cin >> n;
    for (int i = 1; i <= n; ++i) {
        cin >> A[i] >> B[i];
    }
    cout << excrt(n);
    return 0;
}
```



### 倍增求矩阵幂和

```cpp
struct mat {
    int data[N][N] = {};
    int size{};

    int *operator[](int index) {
        return data[index];
    }
} G, e;

mat operator+(const mat &a, const mat &b) {
    mat ret;
    ret.size = a.size;
    for (int i = 1; i <= a.size; ++i) {
        for (int j = 1; j <= a.size; ++j) {
            ret.data[i][j] = a.data[i][j] + b.data[i][j];
        }
    }
    return ret;
}

mat mul(mat &A, mat &B) {
    mat C;
    C.size = A.size;
    for (int i = 1; i <= A.size; i++) {
        for (int k = 1; k <= A.size; k++) {
            for (int j = 1; j <= A.size; j++) {
                C[i][j] = (C[i][j] + A[i][k] * B[k][j]) % mod;
            }
        }
    }
    return C;
}

mat matpow(mat A, int n) {
    mat B;
    B.size = A.size;
    for (int i = 1; i <= A.size; i++) {
        B[i][i] = 1;
    }
    while (n) {
        if (n & 1) B = mul(B, A);
        A = mul(A, A);
        n >>= 1;
    }
    return B;
}

/*倍增法求解A^1 + A^2 + ... + A^n*/
mat pow_sum(const mat &a, int n) {
    if (n == 1) return a;
    mat tmp = pow_sum(a, n / 2);
    mat tt = matpow(a, n / 2);
    mat sum = tmp + mul(tmp, tt);
    ///若n为奇数，n/2 + n/2 = n-1, 因此sum需要加上A^(n)这一项
    if (n & 1) sum = sum + matpow(a, n);
    return sum;
}
```



### 高斯消元

```cpp
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

const int maxn = 10;
const int maxq = 100;

int n;
double a[maxn][maxn + 1]; //系数矩阵
double res[maxn] = {}; //结果数组
double query[maxq];

void solve() {
    int now; // 当前所在行
    double temp; //用于记录消元时的因数
    for (int j = 0; j < n; j++) {
        now = j;
        // 找到当前列最大系数
        for (int i = j; i < n; i++)
            if (fabs(a[i][j]) > fabs(a[now][j]))
                now = i;
        if (now != j)
            //行交换
            for (int i = 0; i < n + 1; i++) {
                double t = a[j][i];
                a[j][i] = a[now][i];
                a[now][i] = t;
            }
        // 消元
        for (int i = j + 1; i < n; i++) {
            temp = a[i][j] / a[j][j];
            for (int k = j; k <= n; k++)
                a[i][k] -= a[j][k] * temp;
        }
    }
    //获得结果
    for (int i = n - 1; i >= 0; i--) {
        if (a[i][i] == 0) {
            res[i] = 0;
            continue;
        }
        res[i] = a[i][n];
        for (int j = i + 1; j < n; j++) {
            res[i] -= a[i][j] * res[j];
        }
        res[i] /= a[i][i];
    }
}

int main() {
    cin >> n;
    if (n == 1) cin >> a[0][0] >> a[0][1];
    else {
        for (int i = 0; i < n; i++) {
            double x, y;
            cin >> x >> y;
            a[i][n - 1] = 1; //常数项系数
            a[i][n] = y;
            for (int j = 0; j < n - 1; j++) {
                a[i][j] = pow(x, n - j - 1);
            }
        }
    }
    int q;
    cin >> q;
    for (int i = 0; i < q; i++) cin >> query[i];
    if (n == 1) {
        double k = a[0][1] / a[0][0];
        for (int i = 0; i < q; i++) {
            double ans = query[i] * k;
            if (ans >= -0.005 && ans < 0) ans = 0;
            cout << fixed << setprecision(2) << ans << '\n';
        }
    } else {
        solve();
        for (int i = 0; i < q; i++) {
            double ans = 0;
            for (int j = 0; j < n; j++) {
                ans += pow(query[i], n - j - 1) * res[j];
            }
            if (ans >= -0.005 && ans < 0) ans = 0;
            cout << fixed << setprecision(2) << ans << '\n';
        }
    }
    return 0;
}

```



### Miller-Rabin & Pollard Rho

```cpp
#include <bits/stdc++.h>

using namespace std;
typedef long long ll;

const int Times = 10;
const int N = 5500;

ll ct;
ll fac[N];

ll gcd(ll a, ll b) {
    return b ? gcd(b, a % b) : a;
}

ll multi(ll a, ll b, ll m) {
    ll ret = 0;
    a %= m;
    while (b) {
        if (b & 1) {
            ret = (ret + a) % m;
        }
        b >>= 1;
        a = (a + a) % m;
    }
    return ret;
}

ll qpow(ll a, ll b, ll m) {
    ll ret = 1;
    a %= m;
    while (b) {
        if (b & 1) {
            ret = ret * a % m;
        }
        b >>= 1;
        a = a * a % m;
    }
    return ret;
}

bool Miller_Rabin(ll n) {
    if (n == 2) return true;
    if (n < 2 || !(n & 1)) return false;
    ll m = n - 1;
    int k = 0;
    while ((m & 1) == 0) {
        k++;
        m >>= 1;
    }
    for (int i = 0; i < Times; ++i) {
        ll a = rand() % (n - 1) + 1;
        ll x = qpow(a, m, n);
        ll y = 0;
        for (int j = 0; j < k; ++j) {
            y = multi(x, x, n);
            if (y == 1 && x != 1 && x != n - 1) return false;
            x = y;
        }
        if (y != 1) return false;
    }
    return true;
}

ll pollard_rho(ll n, ll c) {
    ll i = 1, k = 2;
    ll x = rand() % (n - 1) + 1;
    ll y = x;
    while (true) {
        i++;
        x = (multi(x, x, n) + c) % n;
        ll d = gcd((y - x + n) % n, n);
        if (1 < d && d < n) return d;
        if (y == x) return n;
        if (i == k) {
            y = x;
            k <<= 1;
        }
    }
}

void find(ll n, int c) {
    if (n == 1) return;
    if (Miller_Rabin(n)) {
        fac[ct++] = n;
        return;
    }
    ll p = n;
    ll k = c;
    while (p >= n) p = pollard_rho(p, c--);
    find(p, k);
    find(n / p, k);
}

int main() {
    ll n;
    cin >> n;
    find(n, 120);
    sort(fac, fac + ct);
    //排好序的所有质因子 例如60被拆解为 2 2 3 5
}
```



### 矩阵快速幂

```cpp
//矩阵快速幂
#include <vector>
#include <iostream>

using namespace std;
typedef long long ll;
const int maxn = 250;

struct mat {
    ll data[maxn][maxn] = {};
    int size;

    ll *operator[](int index) {
        return data[index];
    }
};

const ll mod = 1e9 + 7;

mat mul(mat &A, mat &B) {     //矩阵乘法
    mat C;
    C.size = A.size;
    for (int i = 0; i < A.size; i++) {
        for (int k = 0; k < A.size; k++) {
            for (int j = 0; j < A.size; j++) {
                C[i][j] = (C[i][j] + A[i][k] * B[k][j]) % mod;
            }
        }
    }
    return C;
}

mat matpow(mat A, ll n) {       //矩阵快速幂
    mat B;
    B.size = A.size;
    for (int i = 0; i < A.size; i++) {
        B[i][i] = 1;
    }
    while (n) {
        if (n & 1) B = mul(B, A);
        A = mul(A, A);
        n >>= 1;
    }
    return B;
}

int main() {
    int n;
    ll k;
    cin >> n >> k;
    mat A;
    A.size = n;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> A[i][j];
        }
    }
    A = matpow(A, k);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    return 0;
}

```



### 快速判断是否能被某个数整除

　（1）被**2**整除的数的特征：一个整数的末位是偶数（0、2、4、6、8）的数能被2整除。

　（2）被**3**整除的数的特征：一个整数的数字和能被3整除，则这个数能被3整除。

　（3）被**4**整除的数的特征：一**个整数的末尾两位数能被4整除则这个数能被4整**除。可以这样快速判断：最后两位数，要是十位是单数，个位就是2或6，要是十位是双数，个位就是0、4、8。

　（4）被**5**整除的数的特征：一个整数的末位是0或者5的数能被5整除。

　（5）被**6**整除的数的特征：一个整数能被2和3整除，则这个数能被6整除。

　（6）被**7**整除的数的特征：“割减法”。若一个整数的个位数字截去，再从余下的数中，减去个位数的2倍，这样，一次次下去，直到能清楚判断为止，如果差是7的倍数（包括0），则这个数能被7整除。过程为：截尾、倍大、相减、验差。

　　例如，判断133是否7的倍数的过程如下：13－3×2＝7，所以133是7的倍数；又例如判断6139是否7的倍数的过程如下：613－9×2＝595 ， 59－5×2＝49，所以6139是7的倍数，余类推。

　（7）被**8**整除的数的特征：一个整数的未尾三位数能被8整除，则这个数能被8整除。

　（8）被**9**整除的数的特征：一个整数的数字和能被9整除，则这个数能被9整除。

　（9）被**10**整除的数的特征：一个整数的末位是0，则这个数能被10整除。

　（10）被**11**整除的数的特征：“奇偶位差法”。一个整数的奇位数字之和与偶位数字之和的差是11的倍数（包括0），则这个数能被11整除。（隔位和相减）

　　例如，判断491678能不能被11整除的过程如下：奇位数字的和9+6+8=23，偶位数位的和4+1+7=12。23-12=11。因此491678能被11整除。

　（11）被**12**整除的数的特征：一个整数能被3和4整除，则这个数能被12整除。

　（12）被**13**整除的数的特征：若一个整数的个位数字截去，再从余下的数中，加上个位数的4倍，这样，一次次下去，直到能清楚判断为止，如果是13的倍数（包括0），则这个数能被13整除。过程为：截尾、倍大、相加、验差。

　（13）被**17**整除的数的特征：若一个整数的个位数字截去，再从余下的数中，减去个位数的5倍，这样，一次次下去，直到能清楚判断为止，如果差是17的倍数（包括0），则这个数能被17整除。过程为：截尾、倍大、相减、验差。

　（14）被**19**整除的数的特征：若一个整数的个位数字截去，再从余下的数中，加上个位数的2倍，这样，一次次下去，直到能清楚判断为止，如果是19的倍数（包括0），则这个数能被19整除。过程为：截尾、倍大、相加、验差。

　（15）被**7、11、13** 整除的数的共同特征：若一个整数的末3位与末3位以前的数字所组成的数之差（以大减小）能被7、11、13 整除，则这个数能被7、11、13 整除。

　　例如：128114，由于128-114=14，14是7的倍数，所以128114能被7整除。64152，由于152-64=88，88是11的倍数，所以64152能被11整除。94146，由于146-94=52，52是13的倍数，所以94146能被13整除。



### 欧拉降幂

![欧拉降幂](C:\Users\Boboge\Documents\GitHub\template\math\欧拉降幂.png)



### 判断第二类斯特林数奇偶性

```cpp
#include <bits/stdc++.h>

using namespace std;
typedef long long ll;

int solve(ll n) {
    int ret = 0;
    while (n) {
        n >>= 1;
        ret += n;
    }
    return ret;
}

int main() {
    ll n, m;
    cin >> n >> m;
    ll a = n - m;
    ll b = (m - 1) / 2;
    ll ans1 = solve(a + b);
    ll ans2 = solve(a);
    ll ans3 = solve(b);
    if (ans1 == ans2 + ans3) {
        printf("1\n");
    } else {
        printf("0\n");
    }
    return 0;
}
```



### 三分

##### 凸函数

整数

```cpp
//L, R保证在凸函数两端
while (l + 1 < r) {
    int lm = (l + r) >> 1, rm = (lm + r) >> 1;
    if (calc(lm) > calc(rm))
        r = rm;
    else
        l = lm;
}
//答案取 L
```

double

```cpp
while (l + eps < r) {
    double lm = (l + r) / 2, rm = (lm + r) / 2;
    if (calc(lm) > calc(rm))
        r = rm;
    else
        l = lm;
}
//答案取 (L + R) / 2
```

##### 凹函数

只需要将check时的符号互换。



### 线性基

```cpp
//线性基 异或和最大
#define ll long long
const int maxn = 51;
ll d[maxn + 10];

bool add(ll x) {
    for (int i = maxn; i >= 0; --i) {
        if ((x >> i)) {
            if (d[i]) x ^= d[i];
            else {
                d[i] = x;
                return true;
            }
        }
    }
    return false;
}


int main(){
	int n;
	ll ans = 0;
	scanf("%d",&n);
	while(n--){
		ll a;
		scanf("%lld",&a);
		add(a);
	}
	for(int i = maxn;i >= 0; --i){
		if((ans ^ d[i]) > ans){
			ans ^= d[i];
		}
	}
	printf("%lld", ans);
	return 0;
}
```



## 6. String

### 字典树

```cpp
//将某个int的二进制字符串加入字典树用例

int T[N][2];

void add(int x) {
  int now = 0;
  for (int i = 30; i >= 0; --i) {
    int bit = x >> i & 1;
    if (!T[now][bit]) T[now][bit] = ++id;
    now = T[now][bit];
  }
}
```



### KMP

```cpp
//kmp
void GetNextval(char* p, int next[])
{
	int pLen = strlen(p);
	next[0] = -1;
	int k = -1;
	int j = 0;
	while (j < pLen - 1)
	{
		//p[k]表示前缀，p[j]表示后缀  
		if (k == -1 || p[j] == p[k])
		{
			++j;
			++k;
			//较之前next数组求法，改动在下面4行
			if (p[j] != p[k])
				next[j] = k;   //之前只有这一行
			else
				//因为不能出现p[j] = p[ next[j ]]，所以当出现时需要继续递归，k = next[k] = next[next[k]]
				next[j] = next[k];
		}
		else
		{
			k = next[k];
		}
	}

int KmpSearch(char* s, char* p)
{
	int i = 0;
	int j = 0;
	int sLen = strlen(s);
	int pLen = strlen(p);
	int next[pLen+1] = {};
	GetNextval(p,next);
	while (i < sLen && j < pLen)
	{
		//①如果j = -1，或者当前字符匹配成功（即S[i] == P[j]），都令i++，j++    
		if (j == -1 || s[i] == p[j])
		{
			i++;
			j++;
		}
		else
		{
			//②如果j != -1，且当前字符匹配失败（即S[i] != P[j]），则令 i 不变，j = next[j]    
			//next[j]即为j所对应的next值      
			j = next[j];
		}
	}
	if (j == pLen)
		return i - j;
	else
		return -1;
}

```



### ac自动机

```cpp
#include<iostream>
#include<cstring>
#include<queue>
#include<algorithm>

using namespace std;

const int maxn = 1e5;
const int maxs = 1e5;

struct Result {
    int num;
    int pos;
} ans[maxn];//所有单词的出现次数

bool operator<(Result a, Result b) {
    if (a.num != b.num)
        return a.num > b.num;
    else
        return a.pos < b.pos;
}

struct ACMachine {
    struct Tree//字典树
    {
        int fail;//失配指针
        int vis[26];//子节点的位置
        int end;//标记以这个节点结尾的单词编号
    } AC[maxs];//Trie树
    int cnt = 0;//Trie的指针

    inline void Clean(int x) {
        memset(AC[x].vis, 0, sizeof(AC[x].vis));
        AC[x].fail = 0;
        AC[x].end = 0;
    }

    inline void Build(string s, int Num) {
        int l = s.length();
        int now = 0;//字典树的当前指针
        for (int i = 0; i < l; ++i)//构造Trie树
        {
            if (AC[now].vis[s[i] - 'a'] == 0)//Trie树没有这个子节点
            {
                AC[now].vis[s[i] - 'a'] = ++cnt;//构造出来
                Clean(cnt);
            }
            now = AC[now].vis[s[i] - 'a'];//向下构造
        }
        AC[now].end = Num;//标记单词结尾
    }

    void Get_fail()//构造fail指针
    {
        queue<int> Q;//队列
        for (int i = 0; i < 26; ++i)//第二层的fail指针提前处理一下
        {
            if (AC[0].vis[i] != 0) {
                AC[AC[0].vis[i]].fail = 0;//指向根节点
                Q.push(AC[0].vis[i]);//压入队列
            }
        }
        while (!Q.empty())//BFS求fail指针
        {
            int u = Q.front();
            Q.pop();
            for (int i = 0; i < 26; ++i)//枚举所有子节点
            {
                if (AC[u].vis[i] != 0)//存在这个子节点
                {
                    AC[AC[u].vis[i]].fail = AC[AC[u].fail].vis[i];
                    //子节点的fail指针指向当前节点的
                    //fail指针所指向的节点的相同子节点
                    Q.push(AC[u].vis[i]);//压入队列
                } else//不存在这个子节点
                    AC[u].vis[i] = AC[AC[u].fail].vis[i];
                //当前节点的这个子节点指向当前节点fail指针的这个子节点
            }
        }
    }

    int AC_Query(string s)//AC自动机匹配
    {
        int l = s.length();
        int now = 0, ans = 0;
        for (int i = 0; i < l; ++i) {
            now = AC[now].vis[s[i] - 'a'];//向下一层
            for (int t = now; t; t = AC[t].fail)//循环求解
                ans[AC[t].end].num++;
        }
        return ans;
    }
} ACM;

int main() {
    string s[300];
    int n;
    while (true) {
        cin >> n;
        if (n == 0)break;
        ACM.cnt = 0;
        ACM.Clean(0);
        for (int i = 1; i <= n; ++i) {
            cin >> s[i];
            ans[i].num = 0;
            ans[i].pos = i;
            ACM.Build(s[i], i);
        }
        ACM.AC[0].fail = 0;//结束标志
        ACM.Get_fail();//求出失配指针
        cin >> s[0];//文本串
        ACM.AC_Query(s[0]);
        sort(&ans[1], &ans[n + 1]);
        cout << ans[1].num << endl;
        cout << s[ans[1].pos] << endl;
        for (int i = 2; i <= n; ++i) {
            if (ans[i].num == ans[i - 1].num)
                cout << s[ans[i].pos] << endl;
            else
                break;
        }
    }
    return 0;
}
```



### 马拉车

```cpp
//Manacher
string Manacher(string s)
{
    /*改造字符串*/
    string res="$#";
    for(int i=0;i<s.size();++i)
    {
        res+=s[i];
        res+="#";
    }

    /*数组*/
    vector<int> P(res.size(),0);
    int mi=0,right=0;   //mi为最大回文串对应的中心点，right为该回文串能达到的最右端的值
    int maxLen=0,maxPoint=0;    //maxLen为最大回文串的长度，maxPoint为记录中心点

    for(int i=1;i<res.size();++i)
    {
        P[i]=right>i ?min(P[2*mi-i],right-i):1;     //关键句，文中对这句以详细讲解
        
        while(res[i+P[i]]==res[i-P[i]])
            ++P[i];
        
        if(right<i+P[i])    //超过之前的最右端，则改变中心点和对应的最右端
        {
            right=i+P[i];
            mi=i;
        }

        if(maxLen<P[i])     //更新最大回文串的长度，并记下此时的点
        {
            maxLen=P[i];
            maxPoint=i;
        }
    }
    return s.substr((maxPoint-maxLen)/2,maxLen-1);
}
```



### 回文自动机

```cpp
#include<iostream>
#include<cstring>
#include<queue>
#include<algorithm>

using namespace std;

const int maxn = 300;
const int maxs = 2e5;

int ans[maxn];

struct ACMachine {
    struct Tree//字典树
    {
        int fail;//失配指针
        int vis[26];//子节点的位置
        int end;//标记以这个节点结尾的单词编号
    } AC[maxs * 10];//Trie树
    int cnt = 0;//Trie的指针

    inline void init(int x) {
        memset(AC[x].vis, 0, sizeof(AC[x].vis));
        AC[x].fail = 0;
        AC[x].end = 0;
    }

    inline void insert(string s, int id) {
        int l = s.length();
        int now = 0;//字典树的当前指针
        for (int i = 0; i < l; ++i)//构造Trie树
        {
            if (AC[now].vis[s[i] - 'a'] == 0)//Trie树没有这个子节点
            {
                AC[now].vis[s[i] - 'a'] = ++cnt;//构造出来
                init(cnt);
            }
            now = AC[now].vis[s[i] - 'a'];//向下构造
        }
        AC[now].end = id;//标记单词结尾
    }

    void solve()//构造fail指针
    {
        AC[0].fail = 0;//结束标志
        queue<int> Q;//队列
        for (int i = 0; i < 26; ++i)//第二层的fail指针提前处理一下
        {
            if (AC[0].vis[i] != 0) {
                AC[AC[0].vis[i]].fail = 0;//指向根节点
                Q.push(AC[0].vis[i]);//压入队列
            }
        }
        while (!Q.empty())//BFS求fail指针
        {
            int u = Q.front();
            Q.pop();
            for (int i = 0; i < 26; ++i)//枚举所有子节点
            {
                if (AC[u].vis[i] != 0)//存在这个子节点
                {
                    AC[AC[u].vis[i]].fail = AC[AC[u].fail].vis[i];
                    //子节点的fail指针指向当前节点的
                    //fail指针所指向的节点的相同子节点
                    Q.push(AC[u].vis[i]);//压入队列
                } else//不存在这个子节点
                    AC[u].vis[i] = AC[AC[u].fail].vis[i];
                //当前节点的这个子节点指向当
                //前节点fail指针的这个子节点
            }
        }
    }

    int query(string s)//AC自动机匹配
    {
        int l = s.length();
        int now = 0, res = 0;
        for (int i = 0; i < l; ++i) {
            now = AC[now].vis[s[i] - 'a'];//向下一层
            for (int t = now; t; t = AC[t].fail)//循环求解
                ans[AC[t].end]++;
        }
        return res;
    }
} ACM;

string s[maxn];

int main() {
    ios_base::sync_with_stdio(false);
    int n;
    while (true) {
        scanf("%d", &n);
        if (n == 0) break;
        ACM.cnt = 0;
        ACM.init(0);
        for (int i = 1; i <= n; ++i) {
            s[i].resize(maxs);
            scanf("%s", &s[i][0]);
            s[i].resize(strlen(&s[i][0]));
            ans[i] = 0;
            ACM.insert(s[i], i);
        }
        ACM.solve();//求出失配指针
        string str;
        str.resize(1e6);
        scanf("%s", &str[0]);//文本串
        str.resize(strlen(&str[0]));
        ACM.query(str);
        int resm = 0;
        for (int i = 1; i <= n; i++) resm = max(resm, ans[i]);
        printf("%d\n", resm);
        for (int i = 1; i <= n; i++) if (ans[i] == resm) printf("%s\n", s[i].c_str());
    }
    return 0;
}

```



## 7. Others


### 快读
```cpp
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
//int n, m; read(n, m);
```



### Linux对拍脚本

```shell
while true; do
./make>tmp.in #出数据
./tmp<tmp.in>tmp.out #被测程序
./tmp2<tmp.in>tmp2.out #正确（暴力）程序
if diff tmp.out tmp2.out; then #比较两个输出文件
printf AC #结果相同显示AC
else
echo WA #结果不同显示WA，并退出
#cat tmp.out tmp2.out
exit 0
fi #if的结束标志,与C语言相反，0为真
done # while的结束标志

#BY NICK WONG 2014-08-29
#在终端下，进入当前目录，输入"sh ./nick.sh",（其中nick.sh为当前shell脚本名） '#'表示单行注释
#diff在两文件相同时返回空串
```



### Windows对拍脚本

```powershell
@echo off  
:loop  
    rand.exe > data.in
    std.exe < data.in > std.out
    my.exe < data.in > my.out
    fc my.out std.out 
if not errorlevel 1 goto loop  
pause
goto loop
```

