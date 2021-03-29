#include <bits/stdc++.h>
using namespace std;
using ll = long long;
using pii = pair<int, int>;
const int MAXN = 1000005;
const int LOGN = 21;
int logn[MAXN+5] = {0, 0, 1};
// st[i][j].first表示区间[i, i+2^j-1]的最大值
// st[i][j].second表示区间[i, i+2^j-1]的最大值的下标
// 若有多个最大值, 该下标为最左边的最大值的下标
pii st[MAXN][LOGN + 1];
void pre(int n) {
    // 2^x <= n < 2^(x+1) => 2^(x+1) <= 2n < 2^(x+2)
    for (int i = 3; i < MAXN; i++) logn[i] = logn[i/2] + 1;
    for (int j = 1; j <= LOGN; j++) {
        for (int i = 1; i + (1 << j) - 1 <= n; i++) {
            st[i][j] = max(st[i][j - 1], st[i + (1 << (j - 1))][j - 1]);
        }
    }
}
pii query(int l, int r) {
    int s = logn[r-l+1];
    return max(st[l][s], st[r - (1 << s) + 1][s]);
}