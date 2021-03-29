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