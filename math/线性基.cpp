//线性基
#define ll long long
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