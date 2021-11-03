# 数据结构

### 主席树

```cpp
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

