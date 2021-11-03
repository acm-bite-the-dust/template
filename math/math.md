# 数学

### 倍增求矩阵幂和

```cpp
struct mat {
    int data[N][N] = {};
    int size{};

    int *operator[](int index) {
        return data[index];
    }
} G, e;

mat operator+(const mat &a, const mat &b) {
    mat ret;
    ret.size = a.size;
    for (int i = 1; i <= a.size; ++i) {
        for (int j = 1; j <= a.size; ++j) {
            ret.data[i][j] = a.data[i][j] + b.data[i][j];
        }
    }
    return ret;
}

mat mul(mat &A, mat &B) {
    mat C;
    C.size = A.size;
    for (int i = 1; i <= A.size; i++) {
        for (int k = 1; k <= A.size; k++) {
            for (int j = 1; j <= A.size; j++) {
                C[i][j] = (C[i][j] + A[i][k] * B[k][j]) % mod;
            }
        }
    }
    return C;
}

mat matpow(mat A, int n) {
    mat B;
    B.size = A.size;
    for (int i = 1; i <= A.size; i++) {
        B[i][i] = 1;
    }
    while (n) {
        if (n & 1) B = mul(B, A);
        A = mul(A, A);
        n >>= 1;
    }
    return B;
}

/*倍增法求解A^1 + A^2 + ... + A^n*/
mat pow_sum(const mat &a, int n) {
    if (n == 1) return a;
    mat tmp = pow_sum(a, n / 2);
    mat tt = matpow(a, n / 2);
    mat sum = tmp + mul(tmp, tt);
    ///若n为奇数，n/2 + n/2 = n-1, 因此sum需要加上A^(n)这一项
    if (n & 1) sum = sum + matpow(a, n);
    return sum;
}
```



### 二次剩余

```cpp
/*
 * 二次剩余，mod为奇素数时有解
 * 解最多两个，为相反数x, (mod-x)
 * n为0特判
 * 算法返回-1则无解
 * mod为2返回n
 */
#include<bits/stdc++.h>

using namespace std;
typedef long long ll;
ll w;
struct num {
    ll x, y;
};

num mul(num a, num b, ll p) {
    num ans = {0, 0};
    ans.x = ((a.x * b.x % p + a.y * b.y % p * w % p) % p + p) % p;
    ans.y = ((a.x * b.y % p + a.y * b.x % p) % p + p) % p;
    return ans;
}

ll qpow_real(ll a, ll b, ll p) {
    ll ans = 1;
    while (b) {
        if (b & 1)ans = 1ll * ans % p * a % p;
        a = a % p * a % p;
        b >>= 1;
    }
    return ans % p;
}

ll qpow_imag(num a, ll b, ll p) {
    num ans = {1, 0};
    while (b) {
        if (b & 1)ans = mul(ans, a, p);
        a = mul(a, a, p);
        b >>= 1;
    }
    return ans.x % p;
}

ll solve(ll n, ll p) {
    n %= p;
    if (p == 2)return n;
    if (qpow_real(n, (p - 1) / 2, p) == p - 1)return -1;//不存在
    ll a;
    while (1) {
        a = rand() % p;
        w = ((a * a % p - n) % p + p) % p;
        if (qpow_real(w, (p - 1) / 2, p) == p - 1)break;
    }
    num x = {a, 1};
    return qpow_imag(x, (p + 1) / 2, p);
}

int main() {
    srand(time(0));
    int t;
    scanf("%d", &t);
    while (t--) {
        ll n, p;
        scanf("%lld%lld", &n, &p);
        if (!n) {
            printf("0\n");
            continue;
        }
        ll ans1 = solve(n, p), ans2;
        if (ans1 == -1) {
            printf("Hola!\n");
        } else {
            ans2 = p - ans1;
            if (ans1 > ans2)swap(ans1, ans2);
            if (ans1 == ans2)printf("%lld\n", ans1);
            else printf("%lld %lld\n", ans1, ans2);
        }
    }
}
```



### 高斯消元

```cpp
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

const int maxn = 10;
const int maxq = 100;

int n;
double a[maxn][maxn + 1];
double res[maxn] = {};
double query[maxq];

void solve() {
    int now;
    double temp;
    for (int j = 0; j < n; j++) {
        now = j;
        for (int i = j; i < n; i++)
            if (fabs(a[i][j]) > fabs(a[now][j]))
                now = i;
        if (now != j)
            for (int i = 0; i < n + 1; i++) {
                double t = a[j][i];
                a[j][i] = a[now][i];
                a[now][i] = t;
            }
        // œ˚‘™
        for (int i = j + 1; i < n; i++) {
            temp = a[i][j] / a[j][j];
            for (int k = j; k <= n; k++)
                a[i][k] -= a[j][k] * temp;
        }
    }
    for (int i = n - 1; i >= 0; i--) {
        if (a[i][i] == 0) {
            res[i] = 0;
            continue;
        }
        res[i] = a[i][n];
        for (int j = i + 1; j < n; j++) {
            res[i] -= a[i][j] * res[j];
        }
        res[i] /= a[i][i];
    }
}

int main() {
    cin >> n;
    if (n == 1) cin >> a[0][0] >> a[0][1];
    else {
        for (int i = 0; i < n; i++) {
            double x, y;
            cin >> x >> y;
            a[i][n - 1] = 1;
            a[i][n] = y;
            for (int j = 0; j < n - 1; j++) {
                a[i][j] = pow(x, n - j - 1);
            }
        }
    }
    int q;
    cin >> q;
    for (int i = 0; i < q; i++) cin >> query[i];
    if (n == 1) {
        double k = a[0][1] / a[0][0];
        for (int i = 0; i < q; i++) {
            double ans = query[i] * k;
            if (ans >= -0.005 && ans < 0) ans = 0;
            cout << fixed << setprecision(2) << ans << '\n';
        }
    } else {
        solve();
        for (int i = 0; i < q; i++) {
            double ans = 0;
            for (int j = 0; j < n; j++) {
                ans += pow(query[i], n - j - 1) * res[j];
            }
            if (ans >= -0.005 && ans < 0) ans = 0;
            cout << fixed << setprecision(2) << ans << '\n';
        }
    }
    return 0;
}
```



### 矩阵快速幂

```cpp
//矩阵快速幂
#include <vector>
#include <iostream>

using namespace std;
typedef long long ll;
const int maxn = 250;

struct mat {
    ll data[maxn][maxn] = {};
    int size;

    ll *operator[](int index) {
        return data[index];
    }
};

const ll mod = 1e9 + 7;

mat mul(mat &A, mat &B) {     //矩阵乘法
    mat C;
    C.size = A.size;
    for (int i = 0; i < A.size; i++) {
        for (int k = 0; k < A.size; k++) {
            for (int j = 0; j < A.size; j++) {
                C[i][j] = (C[i][j] + A[i][k] * B[k][j]) % mod;
            }
        }
    }
    return C;
}

mat matpow(mat A, ll n) {       //矩阵快速幂
    mat B;
    B.size = A.size;
    for (int i = 0; i < A.size; i++) {
        B[i][i] = 1;
    }
    while (n) {
        if (n & 1) B = mul(B, A);
        A = mul(A, A);
        n >>= 1;
    }
    return B;
}

int main() {
    int n;
    ll k;
    cin >> n >> k;
    mat A;
    A.size = n;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> A[i][j];
        }
    }
    A = matpow(A, k);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }
    return 0;
}
```



### 快速判断能否被一些数整除

1. 被**2**整除的数的特征：一个整数的末位是偶数（0、2、4、6、8）的数能被2整除。

2. 被**3**整除的数的特征：一个整数的数字和能被3整除，则这个数能被3整除。

3. 被**4**整除的数的特征：一**个整数的末尾两位数能被4整除则这个数能被4整**除。可以这样快速判断：最后两位数，要是十位是单数，个位就是2或6，要是十位是双数，个位就是0、4、8。

4. 被**5**整除的数的特征：一个整数的末位是0或者5的数能被5整除。

5. 被**6**整除的数的特征：一个整数能被2和3整除，则这个数能被6整除。

6. 被**7**整除的数的特征：“割减法”。若一个整数的个位数字截去，再从余下的数中，减去个位数的2倍，这样，一次次下去，直到能清楚判断为止，如果差是7的倍数（包括0），则这个数能被7整除。过程为：截尾、倍大、相减、验差。例如，判断133是否7的倍数的过程如下：13－3×2＝7，所以133是7的倍数；又例如判断6139是否7的倍数的过程如下：613－9×2＝595 ， 59－5×2＝49，所以6139是7的倍数，余类推。

7. 被**8**整除的数的特征：一个整数的未尾三位数能被8整除，则这个数能被8整除。

8. 被**9**整除的数的特征：一个整数的数字和能被9整除，则这个数能被9整除。

9. 被**10**整除的数的特征：一个整数的末位是0，则这个数能被10整除。

10. 被**11**整除的数的特征：“奇偶位差法”。一个整数的奇位数字之和与偶位数字之和的差是11的倍数（包括0），则这个数能被11整除。（隔位和相减）。例如，判断491678能不能被11整除的过程如下：奇位数字的和9+6+8=23，偶位数位的和4+1+7=12。23-12=11。因此491678能被11整除。

11. 被**12**整除的数的特征：一个整数能被3和4整除，则这个数能被12整除。

12. 被**13**整除的数的特征：若一个整数的个位数字截去，再从余下的数中，加上个位数的4倍，这样，一次次下去，直到能清楚判断为止，如果是13的倍数（包括0），则这个数能被13整除。过程为：截尾、倍大、相加、验差。

13. 被**17**整除的数的特征：若一个整数的个位数字截去，再从余下的数中，减去个位数的5倍，这样，一次次下去，直到能清楚判断为止，如果差是17的倍数（包括0），则这个数能被17整除。过程为：截尾、倍大、相减、验差。

14. 被**19**整除的数的特征：若一个整数的个位数字截去，再从余下的数中，加上个位数的2倍，这样，一次次下去，直到能清楚判断为止，如果是19的倍数（包括0），则这个数能被19整除。过程为：截尾、倍大、相加、验差。

15. 被**7、11、13** 整除的数的共同特征：若一个整数的末3位与末3位以前的数字所组成的数之差（以大减小）能被7、11、13 整除，则这个数能被7、11、13 整除。例如：128114，由于128-114=14，14是7的倍数，所以128114能被7整除。64152，由于152-64=88，88是11的倍数，所以64152能被11整除。94146，由于146-94=52，52是13的倍数，所以94146能被13整除。



### 欧拉函数

```cpp
//sqrt求单个欧拉函数
ll ph(ll x) {
    ll ret = x, tmp = x;
    for (ll i = 2; i * i <= x; i++) {
        if (tmp % i == 0) {
            ret = ret / i * (i - 1);
            while (tmp % i == 0) tmp /= i;
        }
    }
    if (tmp > 1) ret = ret / tmp * (tmp - 1);
    return ret;
}


//线性筛欧拉函数
int phi[N];
vector<int> prime;
bool isprime[N];
void init() {
    memset(isprime, 1, sizeof(isprime));
    phi[1] = 1, isprime[1] = false;
    for (int i = 2; i < N; ++i) {
        if (isprime[i]) {
            prime.push_back(i);
            phi[i] = i - 1;
        }
        for (int p : prime) {
            if (1ll * i * p >= N) break;
            int now = i * p;
            isprime[now] = false;
            if (i % p == 0) {
                phi[now] = phi[i] * p;
                break;
            } else {
                phi[now] = phi[i] * (p - 1);
            }
        }
    }
}
```



### 欧拉降幂

![欧拉降幂](欧拉降幂.png)



### 判断第二类斯特林数奇偶性

```cpp
#include <bits/stdc++.h>

using namespace std;
typedef long long ll;

int solve(ll n) {
    int ret = 0;
    while (n) {
        n >>= 1;
        ret += n;
    }
    return ret;
}

int main() {
    ll n, m;
    cin >> n >> m;
    ll a = n - m;
    ll b = (m - 1) / 2;
    ll ans1 = solve(a + b);
    ll ans2 = solve(a);
    ll ans3 = solve(b);
    if (ans1 == ans2 + ans3) {
        printf("1\n");
    } else {
        printf("0\n");
    }
    return 0;
}
```



### 三分

#### 凸函数

整数

```cpp
//L, R保证在凸函数两端
while (l + 1 < r) {
    int lm = (l + r) >> 1, rm = (lm + r) >> 1;
    if (calc(lm) > calc(rm))
        r = rm;
    else
        l = lm;
}
//答案取 L
```

double

```cpp
while (l + eps < r) {
    double lm = (l + r) / 2, rm = (lm + r) / 2;
    if (calc(lm) > calc(rm))
        r = rm;
    else
        l = lm;
}
//答案取 (L + R) / 2
```

#### 凹函数

只需要将check时的符号互换。



### 线性基

```cpp
//线性基build
typedef long long ll;
const int maxn = 50;
ll d[maxn + 1];

bool add(ll x) {
    for (int i = maxn; i >= 0; --i) {
        if ((x >> i)) {
            if (d[i]) x ^= d[i];
            else {
                d[i] = x;
                return true;
            }
        }
    }
    return false;
}
```



### exBSGS

```cpp
/* exBSGS算法
 * N ^ x = M (mod P), 其中N, P互质
 * 返回的x为最小解
 */
#include<bits/stdc++.h>

using namespace std;
typedef long long ll;
unordered_map<ll, int> H;

ll gcd(ll a, ll b) {
    if (!b) return a;
    return gcd(b, a % b);
}

ll expow(ll a, ll b, ll mod) {
    ll res = 1;
    while (b) res = ((b & 1) ? res * a % mod : res), a = a * a % mod, b >>= 1;
    return res;
}

ll exgcd(ll &x, ll &y, ll a, ll b) {
    if (!b) {
        x = 1, y = 0;
        return a;
    }
    ll t = exgcd(y, x, b, a % b);
    y -= x * (a / b);
    return t;
}

ll BSGS(ll a, ll b, ll mod, ll q) {
    H.clear();
    ll Q, p = ceil(sqrt(mod)), x, y;
    exgcd(x, y, q, mod), b = (b * x % mod + mod) % mod,
                         Q = expow(a, p, mod), exgcd(x, y, Q, mod), Q = (x % mod + mod) % mod;
    for (ll i = 1, j = 0; j <= p; ++j, i = i * a % mod) if (!H.count(i)) H[i] = j;
    for (ll i = b, j = 0; j <= p; ++j, i = i * Q % mod) if (H[i]) return j * p + H[i];
    return -1;
}

ll exBSGS(ll N, ll M, ll P) {
    ll q = 1;
    ll k = 0, ret;
    if (M == 1) return 0;
    while ((ret = gcd(N, P)) > 1) {
        if (M % ret) return -1;
        ++k, M /= ret, P /= ret, q = q * (N / ret) % P;
        if (q == M) return k;
    }
    return (ret = BSGS(N, M, P, q)) == -1 ? -1 : ret + k;
}

int main() {
    while (true) {
        int N, M, P;
        scanf("%d%d%d", &N, &P, &M);
        if (!N && !M && !P) break;
        N %= P, M %= P;
        int ans = exBSGS(N, M, P);
        if (ans == -1)
            printf("No Solution\n");
        else
            printf("%d\n", ans);
    }
}
```



### exCRT

```cpp
//x % A[i] = B[i] O(nlogn) x为最小解 x + k * lcm都可行
#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
const int N = 1e5 + 10;

ll mul(ll a, ll b, ll mod) {
    ll ret = 0;
    while (b) {
        if (b & 1) ret = (ret + a) % mod;
        a = (a + a) % mod;
        b >>= 1;
    }
    return ret;
}

ll exgcd(ll a, ll b, ll &x, ll &y) {
    ll ret, tmp;
    if (!b) {
        x = 1;
        y = 0;
        return a;
    }
    ret = exgcd(b, a % b, x, y);
    tmp = x;
    x = y;
    y = tmp - a / b * y;
    return ret;
}

ll A[N], B[N];

ll excrt(int n) {
    ll x, y;
    ll M = A[1], ans = B[1];
    for (int i = 2; i <= n; ++i) {
        ll a = M, b = A[i], c = (B[i] - ans % b + b) % b;
        ll g = exgcd(a, b, x, y), bg = b / g;
        if (c % g) return -1;
        x = mul(x, c / g, bg); //可能溢出
        ans += x * M;
        M *= bg;
        ans = (ans % M + M) % M;
    }
    return (ans % M + M) % M;
}

int main() {
    int n;
    cin >> n;
    for (int i = 1; i <= n; ++i) {
        cin >> A[i] >> B[i];
    }
    cout << excrt(n);
    return 0;
}
```



### Miller-Rabin & Pollard-Rho

```cpp
#include <bits/stdc++.h>

using namespace std;
typedef long long ll;

const int Times = 10;
const int N = 5500;

ll ct;
ll fac[N];

ll gcd(ll a, ll b) {
    return b ? gcd(b, a % b) : a;
}

ll multi(ll a, ll b, ll m) {
    ll ret = 0;
    a %= m;
    while (b) {
        if (b & 1) {
            ret = (ret + a) % m;
        }
        b >>= 1;
        a = (a + a) % m;
    }
    return ret;
}

ll qpow(ll a, ll b, ll m) {
    ll ret = 1;
    a %= m;
    while (b) {
        if (b & 1) {
            ret = ret * a % m;
        }
        b >>= 1;
        a = a * a % m;
    }
    return ret;
}

bool Miller_Rabin(ll n) {
    if (n == 2) return true;
    if (n < 2 || !(n & 1)) return false;
    ll m = n - 1;
    int k = 0;
    while ((m & 1) == 0) {
        k++;
        m >>= 1;
    }
    for (int i = 0; i < Times; ++i) {
        ll a = rand() % (n - 1) + 1;
        ll x = qpow(a, m, n);
        ll y = 0;
        for (int j = 0; j < k; ++j) {
            y = multi(x, x, n);
            if (y == 1 && x != 1 && x != n - 1) return false;
            x = y;
        }
        if (y != 1) return false;
    }
    return true;
}

ll pollard_rho(ll n, ll c) {
    ll i = 1, k = 2;
    ll x = rand() % (n - 1) + 1;
    ll y = x;
    while (true) {
        i++;
        x = (multi(x, x, n) + c) % n;
        ll d = gcd((y - x + n) % n, n);
        if (1 < d && d < n) return d;
        if (y == x) return n;
        if (i == k) {
            y = x;
            k <<= 1;
        }
    }
}

void find(ll n, int c) {
    if (n == 1) return;
    if (Miller_Rabin(n)) {
        fac[ct++] = n;
        return;
    }
    ll p = n;
    ll k = c;
    while (p >= n) p = pollard_rho(p, c--);
    find(p, k);
    find(n / p, k);
}

int main() {
    ll n;
    cin >> n;
    find(n, 120);
    sort(fac, fac + ct);
    //排好序的所有质因子 例如60被拆解为 2 2 3 5
}
```

