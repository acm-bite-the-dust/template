#include <iostream>
#include <ctime>
#include <cstdlib>

using namespace std;
typedef long long ll;
const int count = 10; //随机测试次数

ll qpow(ll a, ll n, ll mod) {
    ll ans = 1;
    while (n) {
        if (n & 1) {
            ans = ans * a % mod;
        }
        n >>= 1;
        a = a * a % mod;
    }
    return ans;
}

bool Miller_Rabin(ll n) {
    if (n == 2) return true;
    for (int i = 0; i < count; ++i) {
        int a = rand() % (n - 2) + 2;
        if (qpow(a, n, n) != a)
            return false;
    }
    return true;
}

int main() {
    srand(time(NULL));
    ll n;
    cin >> n;
    if (Miller_Rabin(n))
        cout << "YES" << '\n';
    else
        cout << "NO" << '\n';
    return 0;
}
