bool cycle;
vector<int> topo;
void dfs(int now) {
	vis[now] = 1;
    for(int i : G[now]) {
        if(vis[i] == 0) {
            dfs(i);
        } else if(vis[i] == 1) {
            cycle = true;
        }
    }
    vis[now] = 2;
    topo.push_back(now);
}

reverse(topo.begin(), topo.end());