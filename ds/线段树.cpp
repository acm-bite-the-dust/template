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
