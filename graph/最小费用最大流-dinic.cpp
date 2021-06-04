namespace MCMF
{
const int MAXN = ___, MAXM = ___, INF = 0x7fffffff;
int head[MAXN], cnt = 1;
struct Edge
{
    int to, w, c, next;
} edges[MAXM * 2];
inline void add(int from, int to, int w, int c)
{
    edges[++cnt] = {to, w, c, head[from]};
    head[from] = cnt;
}
inline void addEdge(int from, int to, int w, int c)
{
    add(from, to, w, c);
    add(to, from, 0, -c);
}
int s, t, dis[MAXN], cur[MAXN];
bool inq[MAXN], vis[MAXN];
queue<int> Q;
bool SPFA()
{
    while (!Q.empty())
        Q.pop();
    copy(head, head + MAXN, cur);
    fill(dis, dis + MAXN, INF);
    dis[s] = 0;
    Q.push(s);
    while (!Q.empty())
    {
        int p = Q.front();
        Q.pop();
        inq[p] = 0;
        for (int e = head[p]; e != 0; e = edges[e].next)
        {
            int to = edges[e].to, vol = edges[e].w;
            if (vol > 0 && dis[to] > dis[p] + edges[e].c)
            {
                dis[to] = dis[p] + edges[e].c;
                if (!inq[to])
                {
                    Q.push(to);
                    inq[to] = 1;
                }
            }
        }
    }
    return dis[t] != INF;
}
int dfs(int p = s, int flow = INF)
{
    if (p == t)
        return flow;
    vis[p] = 1;
    int rmn = flow;
    for (int eg = cur[p]; eg && rmn; eg = edges[eg].next)
    {
        cur[p] = eg;
        int to = edges[eg].to, vol = edges[eg].w;
        if (vol > 0 && !vis[to] && dis[to] == dis[p] + edges[eg].c)
        {
            int c = dfs(to, min(vol, rmn));
            rmn -= c;
            edges[eg].w -= c;
            edges[eg ^ 1].w += c;
        }
    }
    vis[p] = 0;
    return flow - rmn;
}
int maxflow, mincost;
inline void run(int s, int t)
{
    MCMF::s = s, MCMF::t = t;
    while (SPFA())
    {
        int flow = dfs();
        maxflow += flow;
        mincost += dis[t] * flow;
    }
}
} // namespace MCMF