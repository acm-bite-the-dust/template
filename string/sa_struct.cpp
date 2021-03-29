#include <bits/stdc++.h>
using namespace std;
const int N = 1000010;
struct SA {
    char s[N]; // 下标从1开始
    int n, chset_siz; // n: 长度      chset_siz: 字符集大小
    int sa[N], rk[N], oldrk[N*2], id[N], px[N], cnt[N], ht[N];
    bool cmp(int x, int y, int w) {
        return oldrk[x] == oldrk[y] && oldrk[x + w] == oldrk[y + w];
    }

    // s[x, n]和s[y, n]的最长公共前缀的长度
    int lcp(int x, int y) {
        if (x == y) return n + 1 - x;
        int rkA = rk[x], rkB = rk[y];
        if (rkA > rkB) swap(rkA, rkB);
        return min(ht[rkA .. rkB]);
    }

    // s[al, ar] < s[bl, br], return -1
    // s[al, ar] == s[bl, br], return 0
    // s[al, ar] > s[bl, br], return 1
    int substr_cmp(int al, int ar, int bl, int br) {
        int lenA = ar - al + 1, lenB = br - bl + 1;
        int _lcp = lcp(al, bl);
        if (lenA == lenB && _lcp >= lenA) return 0;
        if (_lcp > min(lenA, lenB)) return lenA < lenB ? -1 : 1;
        else return rk[al] < rk[bl] ? -1 : 1;
    }

    // 本质不同的子串个数
    int diff_substr_cnt() {
        int res = n * (n + 1) / 2;
        for (int i = 2; i <= n; i++) res -= ht[i];
        return res;
    }

    // 长度 <= k 的本质不同的子串的个数
    int lenK_substr_cnt(int k) {
        int res = min(k, n+1-sa[1]);
        for (int i = 2; i <= n; i++) res += min(k, n+1-sa[i]) - min(ht[i], k);
        return res;
    }
    void init() {
        chset_siz = 300;
        memset(cnt, 0, (chset_siz+4)*sizeof(int));
        scanf("%s", s+1);
        n = strlen(s+1);
        memset(oldrk, 0, (n*2+4)*sizeof(int));
    }

    // 计算sa[], rk[], ht[]
    void calc() {
        int i, m = chset_siz, p, w;
        for (i = 1; i <= n; ++i) ++cnt[rk[i] = s[i]];
        for (i = 1; i <= m; ++i) cnt[i] += cnt[i - 1];
        for (i = n; i >= 1; --i) sa[cnt[rk[i]]--] = i;
        for (w = 1;; w <<= 1, m = p) {
            for (p = 0, i = n; i > n - w; --i) id[++p] = i;
            for (i = 1; i <= n; ++i)
                if (sa[i] > w) id[++p] = sa[i] - w;
            memset(cnt, 0, (m+4)*sizeof(int));
            for (i = 1; i <= n; ++i) ++cnt[px[i] = rk[id[i]]];
            for (i = 1; i <= m; ++i) cnt[i] += cnt[i - 1];
            for (i = n; i >= 1; --i) sa[cnt[px[i]]--] = id[i];
            memcpy(oldrk, rk, sizeof(rk));
            for (p = 0, i = 1; i <= n; ++i)
                rk[sa[i]] = cmp(sa[i], sa[i - 1], w) ? p : ++p;
            if (p == n) break;
        }
        for (i = 1, p = 0; i <= n; ++i) {
            if (p) --p;
            while (s[i + p] == s[sa[rk[i] - 1] + p]) ++p;
            ht[rk[i]] = p;
        }
    }
} sa;

// 长度在区间[L, R]内的本质不同的子串的个数
int main() {
    int T; scanf("%d", &T);
    for (int iT = 1; iT <= T; iT++) {
        sa.init(); sa.calc();
        int l, r; scanf("%d%d", &l, &r);
        int ansL = sa.lenK_substr_cnt(l-1);
        int ansR = sa.lenK_substr_cnt(r);
        printf("Case %d: %d\n", iT, ansR-ansL);
    }
    return 0;
}
