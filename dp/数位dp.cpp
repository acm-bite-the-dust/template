// Problem: P2657 [SCOI2009] windy 数 相邻数位差至少为2的个数，不含前导零，[l,r]个数
// Contest: Luogu
// URL: https://www.luogu.com.cn/problem/P2657
// Memory Limit: 125 MB
// Time Limit: 1000 ms
//
// Powered by CP Editor (https://cpeditor.org)

#include <bits/stdc++.h>

#define pb push_back
#define eb emplace_back
#define size(a) (int)((a).size())

using namespace std;
typedef long long ll;
typedef long double ld;
typedef unsigned long long ull;
const int N = 2e5 + 10;
int dp[10][10][2][2];  // len dig lim st
int num[10];

int calc(int x) {
  if (!x) return 0;
  memset(dp, 0, sizeof(dp));
  int tmp = x;
  int len = 0;
  for (int i = 0; i < 10; ++i) {
    num[i] = tmp % 10;
    tmp /= 10;
    if (num[i] != 0) len = i;
  }
  for (int i = 1; i < num[len]; ++i) dp[len][i][0][1] = 1;
  dp[len][num[len]][1][1] = 1;
  for (int i = len - 1; i >= 0; --i) {
    for (int j = 1; j <= 9; ++j) {
      dp[i][j][0][1] = 1;
    }
    for (int j = 0; j <= 9; ++j) {
      for (int k = 0; k <= 9; ++k) {
        if (abs(j - k) < 2) continue;
        if (k == num[i]) {
          dp[i][k][1][0] += dp[i + 1][j][1][1];
          dp[i][k][1][0] += dp[i + 1][j][1][0];
        } else if (k < num[i]) {
          dp[i][k][0][0] += dp[i + 1][j][1][1];
          dp[i][k][0][0] += dp[i + 1][j][1][0];
        }
        dp[i][k][0][0] += dp[i + 1][j][0][1];
        dp[i][k][0][0] += dp[i + 1][j][0][0];
      }
    }
  }
  int ans = 0;
  for (int i = 0; i <= 9; ++i) {
    for (int j = 0; j <= 1; ++j) {
      for (int k = 0; k <= 1; ++k) {
        ans += dp[0][i][j][k];
      }
    }
  }
  return ans;
}

void solve() {
  int a, b;
  cin >> a >> b;
  int ans = calc(b) - calc(a - 1);
  cout << ans << '\n';
}

int main() {
  ios_base::sync_with_stdio(false);
  int T = 1;
  // cin >> T;
  while (T--) {
    solve();
  }
  return 0;
}