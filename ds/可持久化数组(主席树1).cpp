#include<bits/stdc++.h>

using namespace std;
const int N = 1e6;  // 数据范围
int tot, n, m;
int val[(N << 5) + 10], rt[N + 10], ls[(N << 5) + 10], rs[(N << 5) + 10];
int a[N + 10], len;

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


int build(int l, int r) {  // 建树
    int root = ++tot;
    if (l == r) {
        val[root] = a[l];
        return root;
    }
    int mid = l + r >> 1;
    ls[root] = build(l, mid);
    rs[root] = build(mid + 1, r);
    return root;  // 返回该子树的根节点
}

int update(int k, int pos, int l, int r, int pre) {  // 插入操作
    int dir = ++tot;
    ls[dir] = ls[pre], rs[dir] = rs[pre];
    if (l == r) {
        val[dir] = k;
        return dir;
    }
    int mid = l + r >> 1;
    if (pos <= mid)
        ls[dir] = update(k, pos, l, mid, ls[pre]);
    else
        rs[dir] = update(k, pos, mid + 1, r, rs[pre]);
    return dir;
}

int query(int id, int l, int r, int pos) {  // 查询操作
    if (l == r) return val[id];
    int mid = l + r >> 1;
    if (pos <= mid) {
        return query(ls[id], l, mid, pos);
    } else {
        return query(rs[id], mid + 1, r, pos);
    }
}

int main() {
    read(n, m);
    for (int i = 1; i <= n; ++i) {
        read(a[i]);
    }
    rt[0] = build(1, n);
    int nowv = 1;
    while (m--) {
        int version, op;
        read(version, op);
        if (op == 1) {
            int pos, k;
            read(pos, k);
            rt[nowv++] = update(k, pos, 1, n, rt[version]);
        } else {
            int pos;
            read(pos);
            int ans = query(rt[version], 1, n, pos);
            printf("%d\n", ans);
            rt[nowv++] = rt[version];
        }
    }
    return 0;
}