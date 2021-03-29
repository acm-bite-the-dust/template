/**
 *
 *  T(n) = n * logn
 *  n: sort once      logn: the times of sort
 *
 *  M(n) = s[n] + rk[n] + sa[2n] + oldrk[2n]
 *
 */

#include <bits/stdc++.h>
using namespace std;
const int N = 1000010;

char s[N]; // 下标从1开始
int n, sa[N], rk[N], oldrk[N << 1], id[N], px[N], cnt[N], ht[N];
// px[i] = rk[id[i]]（用于排序的数组所以叫 px）将 rk[id[i]] 存下来，减少不连续内存访问
// 为了防止访问 oldrk[i+w] 导致数组越界，开两倍数组。
// 当然也可以在访问前判断是否越界，但直接开两倍数组方便一些。
// rk[n+1~2n] 应该初始化为0，因为aa < aabb
// 多组的时候要注意，可能上一组rk[1-10]，下一组rk[1-5]，下一组的rk[6-10]需要初始化为0

// s[L1, L1+2x-1]与s[L2, L2+2x-1]的大小
// 可通过比较s[L1, L1+x-1], s[L2, L2+x-1], s[L1+x, L1+2x-1], s[L2+x, L2+x-1]确定
bool cmp(int x, int y, int w) {
    return oldrk[x] == oldrk[y] && oldrk[x + w] == oldrk[y + w];
}

// longest common prefix
// lcp(sa[x], sa[y]) = min(ht[x+1 .. y])
// lcp(x, y) = min(ht[rk[x]+1 .. rk[y]])
// 如果 ht[i] 一直大于某个数，前这么多位就一直没变过; 反之，由于后缀已经排好序了，不可能变了之后变回来
// 举例:
// ht[x-2] = 7, ht[x-1] = 10, ht[x] = 6, ht[x+1] = 11
// 第x-1小后缀与第x-2小后缀的lcp为10, 第x小后缀与第x-1小后缀的lcp为6
// 说明第x小后缀同第x-1小后缀第7位不同, 第x-1小后缀同第x-2小后缀第7位相同
// 所以第x小后缀与第x-2小后缀的lcp为6
// 第x+1小后缀同第x小后缀的lcp为11，两后缀第7位相同
// 所以第x+1小后缀与第x-2小后缀的第7位不同, 两后缀的lcp为6
int lcp(int x, int y) {
    if (x == y) return n + 1 - x;
    int rkA = rk[x], rkB = rk[y];
    if (rkA > rkB) swap(rkA, rkB);
    return min(ht[rkA .. rkB]);
}

// s[al..ar] < s[bl..br], return -1
// s[al..ar] == s[bl..br], return 0
// s[al..ar] > s[bl..br], return 1
int substr_cmp(int al, int ar, int bl, int br) {
    int lenA = ar - al + 1, lenB = br - bl + 1;
    int _lcp = lcp(al, bl);
    if (lenA == lenB && _lcp >= lenA) return 0;
    if (_lcp > min(lenA, lenB)) return lenA < lenB ? -1 : 1;
    else return rk[al] < rk[bl] ? -1 : 1;
}

// 子串就是后缀的前缀，所以可以枚举每个后缀，计算前缀总数，再减掉重复
// 不考虑相同子串出现在不同位置， 长度为n的字符串s有n*(n+1)/2个子串
// 按字典序枚举后缀, 每次添加新后缀的时候统计该后缀的所有前缀中，出现过的前缀的个数
// 即添加第i小后缀使仅考虑s[sa[i]:sa[i]], s[sa[i]:sa[i]+1], s[sa[i]:sa[i]+2]...中出现过的个数
// 1. 添加第1小后缀必不会出现过的前缀
// 2. 添加第i小后缀时, 有ht[i]个前缀曾经出现过，需要减去
//      第i-1小: AaX
//      第  i小: AbY
//    其中A为两后缀的LCP, |A| == ht[i]
//    显然第i小后缀的所有前缀中，长度 <= ht[i]的前缀都在第i-1小后缀中作为前缀出现过
//    由于 a != b, 所以第i-1小的后缀中所有长度大于ht[i]的前缀都没有在第i-1小后缀中作为前缀出现过
//    由于 a < b, 我们是按字典序枚举后缀的，所以第i-1小的后缀中所有长度大于ht[i]的前缀都不会在之前的后缀中作为前缀出现
// 至此, 每个后缀都可以将其前缀划分为这样的两部分: AAAAAABBBBBBB
// 一种是A结尾的前缀，这种前缀只在前一个后缀中作为前缀出现过, 且只出现过一次, 这种前缀的个数是ht[i]
// 一种是B结尾的前缀，这种前缀不会在之前任何一个后缀中作为前缀出现
int diff_substr_cnt() {
    int res = n * (n + 1) / 2;
    for (int i = 2; i <= n; i++) res -= ht[i];
    return res;
}

// 求长度小于等于k的子串个数
int lenK_substr_cnt(int k) {
    int res = min(k, n+1-sa[1]);
    for (int i = 2; i <= n; i++) res += min(k, n+1-sa[i]) - min(ht[i], k);
    return res;
}

// 初始化仅需将cnt[0..|Σ|]置0即可
void init(int m) {
    memset(cnt, 0, (m+4)*sizeof(int));
}

int main() {
    // m是字符集大小
    int i, m = 300, p, w;
    init(m);

    scanf("%s", s + 1);
    n = strlen(s + 1); // n = len(s)
    for (i = 1; i <= n; ++i) ++cnt[rk[i] = s[i]];
    for (i = 1; i <= m; ++i) cnt[i] += cnt[i - 1];
    for (i = n; i >= 1; --i) sa[cnt[rk[i]]--] = i;

    // m=p 就是优化计数排序值域
    // 因为m是字符集大小, p是排序后最大值的最大编号, 即离散化后的集合大小
    for (w = 1;; w <<= 1, m = p) { // w=x表示已知s[i, i+x-1]的sa, rk, 求s[i, i+2x-1]的sa, rk
        // 每次排序用的是基数排序，对基数排序的某个关键字用的是计数排序
        // 因为基数排序所以先第二关键字稳定排序，再第一关键字稳定排序

        // 排序第二关键字
        // id[i]==j, 表示上一趟排序后, 后缀sa[j]按第二关键字排序后排第i位, 可以理解为临时的新sa[]
        for (p = 0, i = n; i > n - w; --i) id[++p] = i; //把s[i+w, i+2w-1]为空串（即第二关键字为无穷小）的位置放前面
        // 把剩下的按排好的顺序放进去
        // 已知后缀sa[i]的大小排在第i, s[sa[i], sa[i]+w-1]是s[sa[i]-w, sa[i]+w-1]的第二关键字
        // sa[i] <= w时，以s[sa[i], sa[i]+w-1]为第二关键字的s[sa[i]-w, sa[i]+w-1]不存在
        for (i = 1; i <= n; ++i)
            if (sa[i] > w) id[++p] = sa[i] - w;

        // 排序第一关键字
        memset(cnt, 0, sizeof(cnt));
        for (i = 1; i <= n; ++i) ++cnt[px[i] = rk[id[i]]];
        for (i = 1; i <= m; ++i) cnt[i] += cnt[i - 1];
        for (i = n; i >= 1; --i) sa[cnt[px[i]]--] = id[i];

        // calc rk[]
        memcpy(oldrk, rk, sizeof(rk)); // 由于计算 rk 的时候原来的 rk 会被覆盖，要先复制一份
        for (p = 0, i = 1; i <= n; ++i)
            // 若两个子串相同，它们对应的 rk 也需要相同，所以要去重
            // 中间过程会出现字串相同, 但最终一定各不相同
            rk[sa[i]] = cmp(sa[i], sa[i - 1], w) ? p : ++p;

        if (p == n) break;
    }

    // ht[i] = lcp(sa[i], sa[i-1])
    // ht[rk[i]] >= ht[rk[i-1]]-1
    for (i = 1, p = 0; i <= n; ++i) {
        if (p) --p;
        while (s[i + p] == s[sa[rk[i] - 1] + p]) ++p;
        ht[rk[i]] = p;
    }

    return 0;
}
