#include <bits/stdc++.h>
using namespace std;
using ll = long long;
const int MAXN = 1000005;
const int MOD = 1000000007;
char s[MAXN], t[MAXN]; // 下标从1开始
int nxt[MAXN]; // s[1..i]的最长公共前后缀的长度, 其实准确的说是pi[]不是nxt[]
int cntTT[MAXN]; // cntTT[i] 表示t[1..i]在t[]中出现的次数
ll cntST[MAXN]; // cntST[|s|+1+i .. |s|+1+|t|] 之和即s的前缀在t中的出现次数

void get_nxt(int len) {
    nxt[1] = 0;
    for (int i = 2; i <= len; i++) {
        int p = nxt[i-1];
        while (p && t[i] != t[p+1]) p = nxt[p];
        nxt[i] = (t[i] == t[p+1] ? p+1 : 0);
    }
}
vector<int> match(int lens, int lent) {
    vector<int> vec;
    int p = 0; // s[?..i-1]的后缀和t已经匹配到p(即已匹配长度)
    for (int i = 1; i <= lens; i++) {
        if (s[i] == t[p+1]) p++;
        else {
            while (p && s[i] != t[p+1]) p = nxt[p];
            if (s[i] == t[p+1]) p++;
        }
        if (p == lent) { // 找到一个匹配
            vec.emplace_back(i-lent+1);
            p = nxt[p];
        }
    }
}
ll prefix_cnt_TT(int len) {
    // todo: make sure that has called get_nxt
    for (int i = 1; i <= len; i++) cntTT[i] = 1;
    for (int i = len; i >= 1; i--) cntTT[nxt[i]] += cntTT[i];
    ll sum = 0;
    for (int i = 1; i <= len; i++) sum += cntTT[i];
    return sum;
}

ll prefix_cnt_ST(int lens, int lent) {
    // todo: make sure that has called get_nxt
    for (int i = 1; i <= lent; i++) cntST[i] = 0;
    int p = 0;
    for (int i = 1; i <= lens; i++) {
        if (s[i] == t[p+1]) {
            p++;
            cntST[p]++;
        } else if (p) {
            p = nxt[p];
            i--;
        }
    }
    for (int i = lent; i >= 1; i--) cntST[nxt[i]] += cntST[i];
    ll sum = 0;
    for (int i = 1; i <= lent; i++) sum += cntST[i];
    return sum;
}

ll pre[MAXN], suf[MAXN];
int main() {
    int T; scanf("%d\n", &T);
    while (T--) {
        scanf("%s%s", s+1, t+1);
        int lens = strlen(s+1);
        int lent = strlen(t+1);
        get_nxt(lent);
        prefix_cnt_ST(lens, lent);
        for (int i = 1; i <= lent; i++) pre[i] = cntST[i];

        reverse(s+1, s+1+lens);
        reverse(t+1, t+1+lent);
        get_nxt(lent);
        prefix_cnt_ST(lens, lent);
        for (int i = 1; i <= lent; i++) suf[i] = cntST[i];

        ll ans = 0;
        for (int i = 1; i < lent; i++) {
            ans += pre[i] * suf[lent - i];
        }
        printf("%lld\n", ans);

    }
}
