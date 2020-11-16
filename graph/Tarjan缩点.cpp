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
        col[u] = num_color;
        ins[u] = false;
        S.pop();
    }
}

/*
 * 染好色后, 同一个color可以缩成一个点, 整张图相当于一张DAG
 * DAG可以拓扑, 搜索啥的, 为所欲为
 */